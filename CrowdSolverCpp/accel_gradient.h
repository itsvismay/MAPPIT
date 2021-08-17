#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

using namespace Eigen;

namespace crowds{
  void accel_gradient(VectorXd& q, int num_agents, int num_points_per_agent, double K_acc, VectorXd& g)
  {
    g.resize(q.size());
    g.setZero();

    double e = 0.0;
    for(int i=0; i<num_agents; i++)
    {
      double mass_i = 1.0;//hard coded for now

      VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
      MatrixXd Q_i = Map<MatrixXd>(q_i.data(), 3, q_i.size()/3).transpose();

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

        Vector3d Fj = K_acc*angle*((v2.cross(z)/v2.squaredNorm()) + (v1.cross(z)/v1.squaredNorm()));
        g(i*3*num_points_per_agent + 3*j+0) = Fj(0);
        g(i*3*num_points_per_agent + 3*j+1) = Fj(1);
        g(i*3*num_points_per_agent + 3*j+2) = Fj(2);
      }

    }
  }
}