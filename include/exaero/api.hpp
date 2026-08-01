#ifndef EXAERO_API_HPP
#define EXAERO_API_HPP

#include <stdint.h>

extern "C" {

int exaero_init(const char* mechanism_name, int32_t total_cells);
int exaero_solve(double* conc_ptr, double* met_ptr, double dt);
void exaero_finalize();

}

#endif // EXAERO_API_HPP
