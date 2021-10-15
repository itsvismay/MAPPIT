#include <Eigen/Dense>

using namespace Eigen;

namespace crowds{

	double reg_energy(VectorXd& q, int num_agents, int num_points_per_agent, VectorXd& K)
	{
		double e = 0.0;
		for(int i=0; i<num_agents; i++)
		{
			VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);

			//size of regular time intervals over the rod;
			double endtime = q_i(q_i.size()-1);
			int segments = num_points_per_agent -1;
			double kt = (endtime/segments);

			double ei =0;
			//For each rod segment, do this:
			for(int s =0; s<num_points_per_agent - 1; s++)
			{
				double x1 = q_i(3*s+0);
				double y1 = q_i(3*s+1);
				double t1 = q_i(3*s+2);
				double x2 = q_i(3*(s+1)+0);
				double y2 = q_i(3*(s+1)+1);
				double t2 = q_i(3*(s+1)+2);

				double dt = t2 - t1 + 1e-6;//add eps in case its 0
				ei += kt/dt;

			}

			e += K(i)*ei;
		}

		return e;
	}
}