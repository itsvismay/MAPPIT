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

#include "../tol_energy.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* variable declarations here */
    
    VectorXd t, UserTols;

    igl::matlab::parse_rhs_double(prhs+0, t);
    igl::matlab::parse_rhs_double(prhs+1, UserTols);
    double K_tol = mxGetScalar(prhs[2]);

    double e_tol = crowds::tol_energy(t, UserTols, K_tol);

    plhs[0] = mxCreateDoubleScalar(e_tol);

}