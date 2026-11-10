# Tau-Argus CSP solver: cell suppression via branch-and-price LP.
#
# The LP back-end is HiGHS (VSCIP target), reached through the JJ* shim in
# src/Jjsolver.c which is compiled against the HiGHS C API. The sources are
# written in C but use C++ I/O (std::cout) and namespaces, so they are
# compiled as C++, exactly as the original NetBeans/MinGW Makefile did with
# g++.

cmake_minimum_required(VERSION 3.16)

# Honor LANGUAGE CXX for legacy .c files without compiler suffix warnings.
if(POLICY CMP0119)
    cmake_policy(SET CMP0119 NEW)
endif()

project(tauargus_csp CXX)

set(CSP_SRC
    src/Cspbridg.c
    src/Cspcard.c
    src/Cspmain.c
    src/Cspnet.c
    src/Jjsolver.c
    src/MT1RC.C
    src/My_time.c
    src/cspaudit.c
    src/cspback.c
    src/cspbranc.c
    src/cspcapa.c
    src/cspcover.c
    src/cspdebug.c
    src/cspgomo.c
    src/cspheur.c
    src/cspprep.c
    src/cspprice.c
    src/cspsep.c
    src/cspsolve.c
)

# All .c files are C++; tell CMake so it picks the C++ compiler.
set_source_files_properties(${CSP_SRC} PROPERTIES LANGUAGE CXX)

# HiGHS version pin. The superbuild exports TAUARGUS_HIGHS_VERSION (read from
# engine/HiGHS-VERSION, the single source of truth). A standalone build from
# the parent checkout reads the same file.
if(NOT DEFINED TAUARGUS_HIGHS_VERSION AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../HiGHS-VERSION")
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/../../HiGHS-VERSION" _TA_HIGHS_TAG)
    string(STRIP "${_TA_HIGHS_TAG}" _TA_HIGHS_TAG)
    string(REGEX REPLACE "^v" "" TAUARGUS_HIGHS_VERSION "${_TA_HIGHS_TAG}")
endif()
find_package(Highs ${TAUARGUS_HIGHS_VERSION} REQUIRED)

option(BUILD_SHARED_LIBS "Build shared native libraries" ON)
add_library(tauargus_csp ${CSP_SRC})
add_library(tauargus::csp ALIAS tauargus_csp)
set(TAUARGUS_MAX_THREADS "" CACHE STRING "Limit HiGHS workers (empty uses its default)")
if(NOT TAUARGUS_MAX_THREADS STREQUAL "")
    if(NOT TAUARGUS_MAX_THREADS MATCHES "^[1-9][0-9]*$")
        message(FATAL_ERROR "TAUARGUS_MAX_THREADS must be a positive integer")
    endif()
    target_compile_definitions(tauargus_csp PUBLIC TAUARGUS_MAX_THREADS=${TAUARGUS_MAX_THREADS})
endif()
if(NOT BUILD_SHARED_LIBS)
    target_compile_definitions(tauargus_csp PUBLIC TAUARGUS_STATIC)
endif()

target_include_directories(tauargus_csp
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        # Callback interfaces (IProgressListener.h, ICallback.h) shared with
        # the hitas module.
        ${CMAKE_CURRENT_SOURCE_DIR}/../hitas/src
)

target_compile_definitions(tauargus_csp PUBLIC VSCIP)

# Not every public function carries an export macro (e.g. TauAuditJj used by
# the Python bindings), so export all symbols from the Windows DLL.
if(WIN32)
    set_target_properties(tauargus_csp PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

target_link_libraries(tauargus_csp PUBLIC highs::highs)

target_compile_options(tauargus_csp PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -fno-strict-aliasing>
)

target_compile_features(tauargus_csp PUBLIC cxx_std_11)

if(WIN32)
    target_compile_definitions(tauargus_csp PRIVATE _CRT_SECURE_NO_WARNINGS MICROSOFT2)
endif()

install(TARGETS tauargus_csp
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Adapter tests (run by the parent repo's dev loop,
# scripts/dev-test.sh). The lifecycle test hooks
# (JJtest_corruptCostShadow, JJtest_shadowCount) are compiled
# into the library only when tests are enabled.
option(TAUARGUS_CSP_BUILD_TESTS "Build csp adapter tests" OFF)
if(TAUARGUS_CSP_BUILD_TESTS)
    enable_testing()
    target_compile_definitions(tauargus_csp PRIVATE TAUARGUS_CSP_UNIT_TEST)
    add_executable(csp_adapter_lifecycle_tests tests/adapter_lifecycle_test.cpp)
    target_link_libraries(csp_adapter_lifecycle_tests PRIVATE tauargus_csp)
    add_executable(csp_adapter_bounds_tests tests/adapter_bounds_test.cpp)
    target_link_libraries(csp_adapter_bounds_tests PRIVATE tauargus_csp)
    add_executable(csp_adapter_feastol_tests tests/adapter_feastol_test.cpp)
    target_link_libraries(csp_adapter_feastol_tests PRIVATE tauargus_csp)
    # Jjsolver.h includes cspmain.h, which needs the hitas callback
    # headers (a PRIVATE include of the library, so not inherited).
    target_include_directories(csp_adapter_lifecycle_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../hitas/src)
    target_include_directories(csp_adapter_bounds_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../hitas/src)
    target_include_directories(csp_adapter_feastol_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../hitas/src)
    add_executable(csp_adapter_infinity_tests tests/adapter_infinity_test.cpp)
    target_link_libraries(csp_adapter_infinity_tests PRIVATE tauargus_csp)
    target_include_directories(csp_adapter_infinity_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../hitas/src)
    add_executable(csp_adapter_infeasible_tests tests/adapter_infeasible_test.cpp)
    target_link_libraries(csp_adapter_infeasible_tests PRIVATE tauargus_csp)
    target_include_directories(csp_adapter_infeasible_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../hitas/src)
    add_test(NAME adapter_lifecycle COMMAND csp_adapter_lifecycle_tests)
    add_test(NAME adapter_bounds COMMAND csp_adapter_bounds_tests)
    add_test(NAME adapter_feastol COMMAND csp_adapter_feastol_tests)
    add_test(NAME adapter_infinity COMMAND csp_adapter_infinity_tests)
    add_test(NAME adapter_infeasible COMMAND csp_adapter_infeasible_tests)
    add_executable(csp_cut_hash_tests tests/cut_hash_test.cpp)
    target_link_libraries(csp_cut_hash_tests PRIVATE tauargus_csp)
    target_include_directories(csp_cut_hash_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../hitas/src)
    add_test(NAME cut_hash COMMAND csp_cut_hash_tests)
    add_executable(csp_audit_fractional_tests tests/audit_fractional_test.cpp)
    target_link_libraries(csp_audit_fractional_tests PRIVATE tauargus_csp)
    add_test(NAME audit_fractional COMMAND csp_audit_fractional_tests)
    add_executable(csp_audit_feastol_tests tests/audit_feastol_test.cpp)
    target_link_libraries(csp_audit_feastol_tests PRIVATE tauargus_csp)
    add_test(NAME audit_feastol COMMAND csp_audit_feastol_tests)
    add_executable(csp_exit_tests tests/csp_exit_test.cpp)
    target_link_libraries(csp_exit_tests PRIVATE tauargus_csp)
    target_include_directories(csp_exit_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../hitas/src)
    add_test(NAME csp_exit COMMAND csp_exit_tests)
    # The forced-failure case must terminate; a run that retries the
    # failed LP solve forever is a failure, so cap the whole test.
    set_tests_properties(csp_exit PROPERTIES TIMEOUT 60)
    add_executable(csp_net_freecol_tests tests/net_freecol_test.cpp)
    target_link_libraries(csp_net_freecol_tests PRIVATE tauargus_csp)
    target_include_directories(csp_net_freecol_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../hitas/src)
    add_test(NAME net_freecol COMMAND csp_net_freecol_tests)
    add_executable(csp_protection_flow_tests tests/protection_flow_test.cpp)
    target_link_libraries(csp_protection_flow_tests PRIVATE tauargus_csp)
    target_include_directories(csp_protection_flow_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../hitas/src)
    add_test(NAME protection_flow COMMAND csp_protection_flow_tests)
    set_tests_properties(protection_flow PROPERTIES TIMEOUT 60)
endif()
