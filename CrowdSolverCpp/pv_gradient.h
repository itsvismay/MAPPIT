#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

using namespace Eigen;

namespace crowds{
  
  void pv_gradient(VectorXd& q, int num_agents, int num_points_per_agent, VectorXd& K_pv, VectorXd& pv, VectorXd& g)
  {
  	g.resize(q.size());
    g.setZero();

	for(int i=0; i<num_agents; i++)
	{

		VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);

		int t_end = q_i.size()-1;
		if(fabs(K_pv(i)) > 1e-6){
			g(i*3*num_points_per_agent + t_end) = K_pv(i)*(q_i(t_end) - pv(i));
		}

		// //For each rod segment, do this:
		// for(int e =0; e<num_points_per_agent - 1; e++)
		// {
		// 	double x1 = q_i(3*e+0);
		// 	double y1 = q_i(3*e+1);
		// 	double t1 = q_i(3*e+2);
		// 	double x2 = q_i(3*(e+1)+0);
		// 	double y2 = q_i(3*(e+1)+1);
		// 	double t2 = q_i(3*(e+1)+2);
		// 	double dt = t2-t1;
		// 	double dx = x2-x1;
		// 	double dy = y2-y1;
		//  double vv = (dt*dt)/(dx*dx + dy*dy);
 			
		// 	double K = K_pv(i)*(vv - pv*pv);

		// 	g(i*3*num_points_per_agent + 3*e+0) += K*((2*dx*dt*dt)/((dx*dx + dy*dy)*(dx*dx + dy*dy))); //dE/dx1
		// 	g(i*3*num_points_per_agent + 3*e+1) += K*((2*dy*dt*dt)/((dx*dx + dy*dy)*(dx*dx + dy*dy))); //dE/dy1
		// 	g(i*3*num_points_per_agent + 3*e+2) += K*(-2*dt/(dx*dx + dy*dy)); //dE/dt1

		// 	g(i*3*num_points_per_agent + 3*(e+1)+0) += K*((-2*dx*dt*dt)/((dx*dx + dy*dy)*(dx*dx + dy*dy))); //dE/dx2
		// 	g(i*3*num_points_per_agent + 3*(e+1)+1) += K*((-2*dx*dt*dt)/((dx*dx + dy*dy)*(dx*dx + dy*dy))); //dE/dy2
		// 	g(i*3*num_points_per_agent + 3*(e+1)+2) += K*(2*dt/(dx*dx + dy*dy)); //dE/dt2

		// }

	}
  }
}