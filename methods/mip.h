#ifndef MIP_H
#define MIP_H

#include "utils.h"
#include <gurobi_c++.h>



ModelResult* solve_mip(ProblemParams const &params, AnnealingState const &state, double timeout);


#endif // MIP_H 
