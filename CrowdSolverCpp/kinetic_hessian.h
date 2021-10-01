#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

using namespace Eigen;
typedef Eigen::Triplet<double> Trip;

namespace crowds{
	void kinetic_hessian(VectorXd& q, int num_agents, int num_points_per_agent, double K_Ke, SparseMatrix<double>& H)
	{
	  std::vector<Trip> Htrips;
	  H.resize(q.size(), q.size());
	  H.setZero();

	  double e = 0.0;
	  for(int i=0; i<num_agents; i++)
	  {
	    double mass_i = 1.0;//hard coded for now

	    VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
	    
	    //E = ((x2 - x1)^2 + (y2 - y1)^2)/(t2 - t1) -> d^2E/dx^2

	    //For each rod segment, do this:
	    for(int e =0; e<num_points_per_agent - 1; e++)
	    {
	      double x1 = q_i(3*e+0);
	      double y1 = q_i(3*e+1);
	      double t1 = q_i(3*e+2);
	      double x2 = q_i(3*(e+1)+0);
	      double y2 = q_i(3*(e+1)+1);
	      double t2 = q_i(3*(e+1)+2);

	      //-----------------
	      //All ddE/dxdy -> 0
	      //-----------------

	      //ddE/dx1dx1, //ddE/dx2dx2
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+0)+0, i*3*num_points_per_agent + 3*(e+0)+0, mass_i/(t2 - t1))); 
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+0, i*3*num_points_per_agent + 3*(e+1)+0, mass_i/(t2 - t1)));

	      //ddE/dy1dy1, //ddE/dy2dy2
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+0)+1, i*3*num_points_per_agent + 3*(e+0)+1, mass_i/(t2 - t1)));  
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+1, i*3*num_points_per_agent + 3*(e+1)+1, mass_i/(t2 - t1)));  

	      //ddE/dx1dx2, ddE/dy1dy2
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+0)+0, i*3*num_points_per_agent + 3*(e+1)+0, -mass_i/(t2 - t1)));  
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+0, i*3*num_points_per_agent + 3*(e+0)+0, -mass_i/(t2 - t1)));  

	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+0)+1, i*3*num_points_per_agent + 3*(e+1)+1, -mass_i/(t2 - t1)));  
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+1, i*3*num_points_per_agent + 3*(e+0)+1, -mass_i/(t2 - t1)));  

	      //-----------------
	      //ddE/dt1dx1, ddE/dx1dt1
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*e+0,     i*3*num_points_per_agent + 3*e+2,-mass_i*(x2 - x1)/((t2-t1)*(t2 - t1)))); 
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*e+2,     i*3*num_points_per_agent + 3*e+0,-mass_i*(x2 - x1)/((t2-t1)*(t2 - t1)))); 


	      //ddE/dt1dx2, ddE/dx2dt1
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+0, i*3*num_points_per_agent + 3*(e+0)+2,	mass_i*(x2 - x1)/((t2-t1)*(t2 - t1)) )     ); 
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+0)+2, i*3*num_points_per_agent + 3*(e+1)+0, mass_i*(x2 - x1)/((t2-t1)*(t2 - t1)) ) ); 

	      //ddE/dt2dx1, ddE/dx1dt2
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+0)+0, i*3*num_points_per_agent + 3*(e+1)+2, mass_i*(x2 - x1)/((t2-t1)*(t2 - t1)) ) ); 
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*(e+0)+0, mass_i*(x2 - x1)/((t2-t1)*(t2 - t1)) ) ); 

	      //ddE/dt2dx2, ddE/dx2dt2
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+0, i*3*num_points_per_agent + 3*(e+1)+2, -mass_i*(x2 - x1)/((t2-t1)*(t2 - t1))) ); 
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*(e+1)+0, -mass_i*(x2 - x1)/((t2-t1)*(t2 - t1))) ); 

	      //-----------------
	      //ddE/dt1dy1, ddE/dy1dt1
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*e+1,     i*3*num_points_per_agent + 3*e+2, -mass_i*(y2 - y1)/((t2-t1)*(t2 - t1)))); 
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*e+2,     i*3*num_points_per_agent + 3*e+1, -mass_i*(y2 - y1)/((t2-t1)*(t2 - t1)))); 

	      //ddE/dt1dy2, ddE/dy2dt1
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+1, i*3*num_points_per_agent + 3*(e+0)+2,mass_i*(y2 - y1)/((t2-t1)*(t2 - t1)) ) ); 
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+0)+2, i*3*num_points_per_agent + 3*(e+1)+1,mass_i*(y2 - y1)/((t2-t1)*(t2 - t1)) ) ); 

	      //ddE/dt2dy1, ddE/dy1dt2
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+0)+1, i*3*num_points_per_agent + 3*(e+1)+2, -mass_i*(y2 - y1)/((t2-t1)*(t2 - t1))) ); 
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*(e+0)+1, -mass_i*(y2 - y1)/((t2-t1)*(t2 - t1))) ); 

	      //ddE/dt2dy2, ddE/dy2dt2
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+1, i*3*num_points_per_agent + 3*(e+1)+2, mass_i*(y2 - y1)/((t2-t1)*(t2 - t1))) ); 
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*(e+1)+1, mass_i*(y2 - y1)/((t2-t1)*(t2 - t1))) ); 

	      //-----------------
	      //ddE/dt1dt1
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*e+2,     i*3*num_points_per_agent + 3*e+2, mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1))/((t2 - t1)*(t2 - t1)*(t2 - t1)))); 
	      
	      //ddE/dt1dt2, ddE/dt1dt2 
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+0)+2, i*3*num_points_per_agent + 3*(e+1)+2, -mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1))/((t2 - t1)*(t2 - t1)*(t2 - t1))) ); 
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*(e+0)+2, -mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1))/((t2 - t1)*(t2 - t1)*(t2 - t1))) ); 

	      //ddE/dt2dt2
	      Htrips.push_back(Trip(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*(e+1)+2, mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1))/((t2 - t1)*(t2 - t1)*(t2 - t1))) ); 
	    }

	  }

	  H.setFromTriplets(Htrips.begin(), Htrips.end());
	  H = K_Ke*H;
	}
}