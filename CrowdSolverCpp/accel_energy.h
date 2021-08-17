#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

using namespace Eigen;

namespace crowds{
  
  double accel_energy(VectorXd& q, int num_agents, int num_points_per_agent, double K_acc)
  {
    double e = 0.0;
    for(int i=0; i<num_agents; i++)
    {
      double mass_i = 1.0;//hard coded for now

      VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
      MatrixXd Q_i = Map<MatrixXd>(q_i.data(), 3, q_i.size()/3).transpose();


      // VectorXd dx = Q_i.col(0).tail(num_points_per_agent -1) - Q_i.col(0).head(num_points_per_agent - 1) ;
      // VectorXd dy = Q_i.col(1).tail(num_points_per_agent -1) - Q_i.col(1).head(num_points_per_agent - 1) ;
      // VectorXd dt = Q_i.col(2).tail(num_points_per_agent -1) - Q_i.col(2).head(num_points_per_agent - 1) ;
      
      double e_i = 0;
      for(int j=1; j<Q_i.cols() -1; j++){
        //for each node thats not the first node (0) and last node
        Vector3d v1 = Q_i.row(j) - Q_i.row(j-1);
        Vector3d v2 = Q_i.row(j+1) - Q_i.row(j);
        Vector3d v1xv2 = (v1.cross(v2)).transpose();
        double v1dv2 = v1.dot(v2);
        Vector3d z = v1xv2/v1xv2.norm();
        double X = v1.norm()*v2.norm() + v1dv2;
        double Y = v1xv2.dot(z);
        double angle = 2*atan2(Y, X);
        e_i += 0.5*K_acc*(angle - 0)*(angle - 0);

      }
      e += e_i;
    }
    return e;
  }
}