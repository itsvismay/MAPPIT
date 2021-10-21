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

#include "../pv_energy.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* variable declarations here */
    
    VectorXd q, K, mass;

    igl::matlab::parse_rhs_double(prhs+0, q);
    igl::matlab::parse_rhs_double(prhs+1, K);
    double pv = mxGetScalar(prhs[2]);
    int num_agents = mxGetScalar(prhs[3]);
    int num_points_per_agent = mxGetScalar(prhs[4]);
    

    double e = crowds::pv_energy(q, num_agents, num_points_per_agent, K, pv);

    plhs[0] = mxCreateDoubleScalar(e);

}