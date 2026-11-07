#include <cstdlib>
#include "Cspdefns.h"
#include "cspsep.h"

/* The cut-pool globals live in the csp library. On Windows that is a DLL
 * whose exported symbols are functions, so the test seeds and reads the
 * pool through these accessors instead of the data symbols. */
void JJtest_seedCutPool(int n, CONSTRAINT **pool, int max_pool, double zero);
int JJtest_cutPoolRows(void);

int main() {
    VARIABLE variables[2] = {};
    variables[0].index = 17;
    variables[1].index = 23;
    VARIABLE *stack[] = {variables, variables + 1};
    double coefficients[] = {1, 2};
    CONSTRAINT original = {};
    original.type = CAPACITY;
    original.card = 2;
    original.stack = stack;
    original.coef = coefficients;
    original.rhs = 1;
    CONSTRAINT *pool[8] = {};
    JJtest_seedCutPool(0, pool, 8, 1e-9);
    std::srand(123);
    const int expected = std::rand();
    std::srand(123);
    if (new_row(&original) != 1) return 1;
    if (std::rand() != expected) return 2;

    CONSTRAINT duplicate = original;
    if (new_row(&duplicate) != 0 || JJtest_cutPoolRows() != 1) return 3;
    VARIABLE *reversed[] = {variables + 1, variables};
    double reversed_coefficients[] = {2, 1};
    duplicate.stack = reversed;
    duplicate.coef = reversed_coefficients;
    if (new_row(&duplicate) != 0 || JJtest_cutPoolRows() != 1) return 4;
    double scaled[] = {2, 4};
    duplicate.stack = stack;
    duplicate.coef = scaled;
    duplicate.rhs = 2;
    if (new_row(&duplicate) != 0 || JJtest_cutPoolRows() != 1) return 5;
    double different[] = {1, 3};
    CONSTRAINT other = original;
    other.coef = different;
    if (new_row(&other) != 1 || JJtest_cutPoolRows() != 2) return 6;
    return 0;
}
