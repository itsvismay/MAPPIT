#include <Eigen/Dense>

using namespace Eigen;

namespace crowds{

	double tol_energy(VectorXd& t, VectorXd& UserTols, double K_tol)
	{
	  return 0.5*K_tol*(t - UserTols).transpose()*(t - UserTols);
	}
}