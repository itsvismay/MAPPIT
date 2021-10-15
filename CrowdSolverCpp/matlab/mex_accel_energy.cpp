//matlab hooks from libigl (taken directly from gptoolbox)
#include <igl/C_STR.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include <igl/list_to_matrix.h>

#include <mex.h>

#include "../accel_energy.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* variable declarations here */
    
    VectorXd q, K;

    igl::matlab::parse_rhs_double(prhs+0, q);
    igl::matlab::parse_rhs_double(prhs+1, K);;
    int num_agents = mxGetScalar(prhs[2]);
    int num_points_per_agent = mxGetScalar(prhs[3]);

    double e_ke = crowds::accel_energy(q, num_agents, num_points_per_agent, K);

    plhs[0] = mxCreateDoubleScalar(e_ke);

}