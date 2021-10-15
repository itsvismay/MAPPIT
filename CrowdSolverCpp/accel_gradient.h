#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#define PI 3.14159265
using namespace Eigen;

namespace crowds{

  // lots of div by 0 causing nans. 
  // if div by 0, then sign should be 0
  double sign(double x1, double x2){
    double s = (x1 - x2)/fabs(x1-x2);
    if (s != s)
      s = 0;

    return s;
  }

  void accel_gradient(VectorXd& q, int num_agents, int num_points_per_agent, VectorXd& K_acc, VectorXd& g)
  {
    g.resize(q.size());
    g.setZero();
    VectorXd g_test = g;

    MatrixXd Q(3*num_points_per_agent, num_agents);
   

    double e = 0.0;
    for(int i=0; i<num_agents; i++)
    {
      double mass_i = 1.0;//hard coded for now

      VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
      MatrixXd Q_i = Map<MatrixXd>(q_i.data(), 3, q_i.size()/3).transpose();


      int j=1;
      for(j=1; j<Q_i.rows() -1; j++)
      {
        Vector3d v1 = Q_i.row(j) - Q_i.row(j-1);
        Vector3d v2 = Q_i.row(j+1) - Q_i.row(j);
        double K = K_acc(i)*(v1.norm() + v2.norm());
        Vector3d v1xv2 = (v1.cross(v2)).transpose();
        double v1dv2 = v1.dot(v2);
        double v1xv2norm = v1xv2.norm();

        // //----Dave's Energy-----
        double eps = 1e-5;
        if(v1dv2<=eps){
          v1dv2 = eps;

        }
       
        v1xv2 = (v1.cross(v2)).transpose();
        v1dv2 = v1.dot(v2);
        v1xv2norm = v1xv2.norm();

        double X = v1dv2;
        double Y = v1xv2norm;

        double dE = 0;
        double datan = 1.0/(1.0 + (Y/X)*(Y/X));
        if(v1xv2.norm()<=eps)
        {
          dE = 0;
          v1xv2norm = eps;
        }
        else
        {
          dE = K*atan2(Y, X);
        }
        double lolo = v1dv2*v1dv2;

        //------------------Analytical gradient
        //dEdv1
        Matrix3d v2CrossMatrix;
        v2CrossMatrix<< 0,  v2(2), -v2(1),
                       -v2(2),   0,  v2(0),
                        v2(1), -v2(0),   0; 
        Vector3d lodhi1 = v1dv2*(1.0/v1xv2norm) *(v1xv2.transpose()*v2CrossMatrix);
        Vector3d hidlo1 = v2*v1xv2norm;

   
        //dEdv2
        Matrix3d v1CrossMatrix;
        v1CrossMatrix<< 0,  -v1(2), v1(1),
                         v1(2),   0,  -v1(0),
                        -v1(1), v1(0),   0; 
        Vector3d lodhi2 = v1dv2*(1.0/v1xv2norm) *(v1xv2.transpose()*v1CrossMatrix);
        Vector3d hidlo2 = v1*v1xv2norm;

        Vector3d Gjprev = -dE*datan*(lodhi1 - hidlo1)/lolo;
        Vector3d Gjnext = dE*datan*(lodhi2 - hidlo2)/lolo;
        Vector3d Gj =  - Gjprev - Gjnext;
        
        if(Gj(0) != Gj(0))
        {
          std::cout<<"NANS"<<std::endl;
          std::cout<<v1xv2norm<<std::endl;
          std::cout<<v1dv2<<std::endl;
          std::cout<<v1.transpose()<<std::endl;
          std::cout<<v2.transpose()<<std::endl;
          exit(0);
        }
        double e_i = 0.5*K*atan2(Y, X)*atan2(Y, X);
        // std::cout<<v1.norm()<<", "<<v2.norm()<<v1dv2<<", "<<v1xv2norm<<", "<< atan2(Y,X)*180.0/PI <<", "<<e_i<<std::endl;

        g(i*3*num_points_per_agent + 3*(j-1)+0) += Gjprev(0);
        g(i*3*num_points_per_agent + 3*(j-1)+1) += Gjprev(1);
        // g(i*3*num_points_per_agent + 3*(j-1)+2) += Gjprev(2);

        g(i*3*num_points_per_agent + 3*j+0)     += Gj(0);
        g(i*3*num_points_per_agent + 3*j+1)     += Gj(1);
        // g(i*3*num_points_per_agent + 3*j+2)     += Gj(2);

        g(i*3*num_points_per_agent + 3*(j+1)+0) += Gjnext(0);
        g(i*3*num_points_per_agent + 3*(j+1)+1) += Gjnext(1);
        // g(i*3*num_points_per_agent + 3*(j+1)+2) += Gjnext(2);

        
      }
      
    }

  }
}