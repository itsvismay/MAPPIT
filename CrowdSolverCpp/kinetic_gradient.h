#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

using namespace Eigen;

namespace crowds{
  void kinetic_gradient(VectorXd& q, int num_agents, int num_points_per_agent, VectorXd& K_Ke, VectorXd& mass, VectorXd& g)
  {
    g.resize(q.size());
    g.setZero();

    double energy = 0.0;
    for(int i=0; i<num_agents; i++)
    {
      double mass_i = mass(i)*K_Ke(i);

      VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
      
      //E = ((x2 - x1)^2 + (y2 - y1)^2)/(t2 - t1) -> dE/dx

      //For each rod segment, do this:
      for(int e =0; e<num_points_per_agent - 1; e++)
      {
        double x1 = q_i(3*e+0);
        double y1 = q_i(3*e+1);
        double t1 = q_i(3*e+2);
        double x2 = q_i(3*(e+1)+0);
        double y2 = q_i(3*(e+1)+1);
        double t2 = q_i(3*(e+1)+2);

        g(i*3*num_points_per_agent + 3*e+0) += -mass_i*(x2 - x1)/(t2 - t1); //dE/dx1
        g(i*3*num_points_per_agent + 3*e+1) += -mass_i*(y2 - y1)/(t2 - t1); //dE/dy1
        g(i*3*num_points_per_agent + 3*e+2) += 0.5*mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) )/((t2 - t1)*(t2 - t1)); //dE/dt1

        g(i*3*num_points_per_agent + 3*(e+1)+0) += mass_i*(x2 - x1)/(t2 - t1); //dE/dx2
        g(i*3*num_points_per_agent + 3*(e+1)+1) += mass_i*(y2 - y1)/(t2 - t1); //dE/dy2
        g(i*3*num_points_per_agent + 3*(e+1)+2) += -0.5*mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) )/((t2 - t1)*(t2 - t1)); //dE/dt2
      }

    }
  }
}