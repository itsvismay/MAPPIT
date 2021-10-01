#include <Eigen/Dense>

using namespace Eigen;

namespace crowds
{
	void reg_gradient(VectorXd& t, VectorXd& UserTols, double K_tol, VectorXd& g)
	{
	  g.resize(t.size());
	  g.setZero();
	  g = K_tol*(t - UserTols);
	}
}