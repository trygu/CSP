# tauargus — csp

The **csp** module is the linear-programming backend of the τ-ARGUS
optimisation step. It implements the classical cell-suppression /
iterative-rounding LP (branch-and-price over the response table) and, since
the native rewrite, the **audit** (intervalle) routine as well.

The LP is solved by **HiGHS**, reached through the `JJ*` shim in
`src/Jjsolver.c` which compiles directly against the HiGHS C API. The sources
are written in C but use C++ I/O and namespaces, so they are compiled as C++.

## What it provides

- **CSP solver** — `Cspmain` / `cspsolve` and friends: the branch-and-price
  LP that finds the protected cell pattern for a response table.
- **Audit (intervalle)** — `cspaudit.c` exposes the native port of the legacy
  `intervalle.exe` audit program as `TauAuditJj` / `TauAuditJjFree`
  (`src/cspaudit.h`). Given a `.jj` table it returns, per cell, the
  audit interval `[min, max]` and whether the cell is **unsafe** (an
  under-protected / primary-unsafe value). This is the engine behind the
  `+AR` audit columns in the Python layer.

The old NetBeans/MinGW build is gone; CMake is the only build system.

## Build

Part of the superbuild in the top-level repo:

```bash
cmake -S native -B native/build \
      -DHighs_DIR=$(brew --prefix)/lib/cmake/Highs
cmake --build native/build
```

Standalone (needs HiGHS):

```bash
cmake -S native/csp -B native/csp/build \
      -DHighs_DIR=$(brew --prefix)/lib/cmake/Highs
cmake --build native/csp/build
```

Produces `libtauargus_csp.dylib` (shared) and links `highs::highs`.

## License

© M.D. Montes de Oca & J.J. Salazar. Distributed under the EUPL v1.1 /
EUPL-1.2. See the top-level repository for the full notice.
