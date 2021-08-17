#include <Eigen/Dense>
using namespace Eigen;

namespace crowds
{
	void tol_hessian(VectorXd& t, VectorXd& UserTols, double K_tol, MatrixXd& H)
	{
	  H.resize(t.size(), t.size());
	  H.setZero();

	  for(int i=0; i<t.size(); i++){
	    H(i,i) = K_tol;
	  }
	}
}