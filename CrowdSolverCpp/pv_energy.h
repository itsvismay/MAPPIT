#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

using namespace Eigen;

namespace crowds{
  
  double pv_energy(VectorXd& q, int num_agents, int num_points_per_agent, VectorXd& K_pv, VectorXd& pv)
  {
	double energy = 0.0;
	for(int i=0; i<num_agents; i++)
	{

		VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);

		// option 1:
		// in this version, pv is the desired end time
		// just get end times to match
		double energy_i = K_pv(i)*0.5*(q_i(q_i.size()-1) - pv(i))*(q_i(q_i.size()-1) - pv(i));
		// std::cout<<"q_i: "<<q_i(q_i.size()-1)<<", "<<pv(i)<<", "<<K_pv(i)<<std::endl;
		if(fabs(K_pv(i))< 1e-6){
			energy_i = 0;
		}

		// option 2, 
		// matching the slope of each rod segment
		//For each rod segment, do this:
		// double energy_i = 0;
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
		// 	double vv = (dt*dt)/(dx*dx + dy*dy);

		// 	energy_i += 0.5*(vv - pv*pv)*(vv - pv*pv);
		// }

		energy += energy_i;
	}
	// std::cout<<"mex: "<<energy<<std::endl;
	return energy;
  }
}