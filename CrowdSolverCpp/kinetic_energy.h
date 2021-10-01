#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

using namespace Eigen;

namespace crowds{
  
  double kinetic_energy(VectorXd& q, int num_agents, int num_points_per_agent, double K_Ke)
  {
    double e = 0.0;
    for(int i=0; i<num_agents; i++)
    {
      double mass_i = 1.0;//hard coded for now

      VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
      MatrixXd Q_i = Map<MatrixXd>(q_i.data(), 3, q_i.size()/3).transpose();
      
      VectorXd dx = Q_i.col(0).tail(num_points_per_agent -1) - Q_i.col(0).head(num_points_per_agent - 1) ;
      VectorXd dy = Q_i.col(1).tail(num_points_per_agent -1) - Q_i.col(1).head(num_points_per_agent - 1) ;
      VectorXd dt = Q_i.col(2).tail(num_points_per_agent -1) - Q_i.col(2).head(num_points_per_agent - 1) ;
      

      ArrayXd dx2dy2_dt = (dx.array()*dx.array() + dy.array()*dy.array())/dt.array();

      double e_i = 0.5*mass_i*dx2dy2_dt.matrix().sum();//e = 0.5*(dx^2 + dy^2)/dt
      e += K_Ke*e_i;
    }
    return e;
  }
}