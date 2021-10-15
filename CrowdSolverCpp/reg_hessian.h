#include <Eigen/Dense>
using namespace Eigen;

typedef Eigen::Triplet<double> Trip;

namespace crowds
{
	void reg_hessian(VectorXd& q, int num_agents, int num_points_per_agent, VectorXd& K, SparseMatrix<double>& H)
	{
		std::vector<Trip> Htrips;
		H.resize(q.size(), q.size());
		H.setZero();

		for(int i=0; i<num_agents; i++)
		{
			VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);

			//size of regular time intervals over the rod;
			double endtime = q_i(q_i.size()-1);
			int segments = num_points_per_agent -1;
			double kt = (endtime/segments);

			//For each rod segment, do this:
			int s = 0;
			for(s =0; s<num_points_per_agent - 1; s++)
			{
				double x1 = q_i(3*s+0);
				double y1 = q_i(3*s+1);
				double t1 = q_i(3*s+2);
				double x2 = q_i(3*(s+1)+0);
				double y2 = q_i(3*(s+1)+1);
				double t2 = q_i(3*(s+1)+2);

				double dt = t2 - t1 + 1e-6;//add eps in case its 0

				Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(s+0)+2,
									  i*3*num_points_per_agent + 3*(s+0)+2, 
									  K(i)*2*kt/(dt*dt*dt))); 

				Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(s+0)+2,
									  i*3*num_points_per_agent + 3*(s+1)+2, 
									  -K(i)*2*kt/(dt*dt*dt))); 

				Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(s+1)+2,
									  i*3*num_points_per_agent + 3*(s+0)+2, 
									  -K(i)*2*kt/(dt*dt*dt)));

				Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(s+1)+2,
									  i*3*num_points_per_agent + 3*(s+1)+2, 
									  K(i)*2*kt/(dt*dt*dt))); 
				

			}

			//Technically the last block should be 0 since 'kt' is not exactly a constant
			//g(i*3*num_points_per_agent + 3*(s+1)+2) = 0;
		}
		H.setFromTriplets(Htrips.begin(), Htrips.end());
	}
}