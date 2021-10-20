#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
using namespace Eigen;

#define PI 3.14159265

namespace crowds{
  
  double accel_energy(VectorXd& q, int num_agents, int num_points_per_agent, VectorXd& K_acc)
  {
    double e = 0.0;
    for(int i=0; i<num_agents; i++)
    {
      VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
      MatrixXd Q_i = Map<MatrixXd>(q_i.data(), 3, q_i.size()/3).transpose();

      double e_i = 0;

      for(int j=1; j<Q_i.rows() -1; j++){
        //for each node thats not the first node (0) and last node
        Vector3d v1 = Q_i.row(j) - Q_i.row(j-1);
        Vector3d v2 = Q_i.row(j+1) - Q_i.row(j);
        double K = K_acc(i)*(v1.norm() + v2.norm());
        Vector3d v1xv2 = v1.cross(v2);
        double v1dv2 = v1.dot(v2);
        double v1xv2norm = v1xv2.norm();
      
        // //----Dave's Energy-----
        double eps = 1e-5;
        if(v1dv2<=eps){
          v1dv2 = eps;

        }
        double X = v1dv2;
        double Y = v1xv2norm;
        e_i += 0.5*K*atan2(Y, X)*atan2(Y, X);

        
      }
      e += e_i;
    }

    return e;
  }
}