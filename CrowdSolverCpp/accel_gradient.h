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

  void accel_gradient(VectorXd& q, int num_agents, int num_points_per_agent, double K_acc, VectorXd& g)
  {
    g.resize(q.size());
    g.setZero();
    VectorXd g_test = g;

    double e = 0.0;
    for(int i=0; i<num_agents; i++)
    {
      double mass_i = 1.0;//hard coded for now

      VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
      MatrixXd Q_i = Map<MatrixXd>(q_i.data(), 3, q_i.size()/3).transpose();
      std::cout<<"Q"<<std::endl;
      std::cout<<Q_i<<std::endl<<std::endl;

      for(int j=1; j<Q_i.rows() -1; j++){
        // //for each node thats not the first node (0) and last node
        // Vector3d v1 = Q_i.row(j) - Q_i.row(j-1);
        // Vector3d v2 = Q_i.row(j+1) - Q_i.row(j);
        // Vector3d v1xv2 = (v1.cross(v2)).transpose();
        // double v1dv2 = v1.dot(v2);
        // double v1xv2norm = v1xv2.norm();
        // double eps = 1e-3;
        // if(v1xv2.norm()<=eps){
        //   v1xv2norm = eps;
        // }
        // Vector3d z = v1xv2/v1xv2norm;
        // double X = v1.norm()*v2.norm() + v1dv2;
        // double Y = v1xv2.dot(z);
        // double angle = 2*atan2(Y, X);

        // Vector3d Gjprev = -K_acc*(angle - 0) * (v1.cross(z)/v1.squaredNorm());
        // Vector3d Gjnext = -K_acc*(angle - 0) * (v2.cross(z)/v2.squaredNorm());

        // //Vector3d Gjprev = angle*(v1.cross(z)/v1.squaredNorm());
        // //Vector3d Gjnext = angle*(v2.cross(z)/v2.squaredNorm());
        // Vector3d Gj =  - Gjprev - Gjnext;

        // g(i*3*num_points_per_agent + 3*(j-1)+0) = Gjprev(0);
        // g(i*3*num_points_per_agent + 3*(j-1)+1) = Gjprev(1);
        // g(i*3*num_points_per_agent + 3*(j-1)+2) = Gjprev(2);

        // g(i*3*num_points_per_agent + 3*j+0) = Gj(0);
        // g(i*3*num_points_per_agent + 3*j+1) = Gj(1);
        // g(i*3*num_points_per_agent + 3*j+2) = Gj(2);

        // g(i*3*num_points_per_agent + 3*(j+1)+0) = Gjnext(0);
        // g(i*3*num_points_per_agent + 3*(j+1)+1) = Gjnext(1);
        // g(i*3*num_points_per_agent + 3*(j+1)+2) = Gjnext(2);

        Vector3d v1 = Q_i.row(j-1) - Q_i.row(j);
        Vector3d v2 = Q_i.row(j+1) - Q_i.row(j);
        double angle = acos(v1.dot(v2)/(v1.norm()*v2.norm()));
        double e_i = 0.5*(angle - PI)*(angle - PI);
        
        double q1 = Q_i.row(j-1)[0];
        double q2 = Q_i.row(j-1)[1];
        double q3 = Q_i.row(j-1)[2];
        double q4 = Q_i.row(j)[0];
        double q5 = Q_i.row(j)[1];
        double q6 = Q_i.row(j)[2];
        double q7 = Q_i.row(j+1)[0];
        double q8 = Q_i.row(j+1)[1];
        double q9 = Q_i.row(j+1)[2];
        double eps = 0;

        /*
        E = 0.5*K*(A - PI)^2
        dE/dx = dE/dA * dA/dv * dv/dx


        */

        // double lo = v1.norm()*v2.norm();
        // double hi = v1.dot(v2);
        // double K = K_acc * (angle-PI) * (-1/sqrt(1 - (hi/lo)*(hi/lo)));
        // //lo'hi - hi'lo /lo*lo
        // Vector3d s1 = ( lo*v1 - hi*v2.norm()*2*v1 )/(lo*lo);
        // Vector3d s3 = ( lo*v2 - hi*v1.norm()*2*v2 ) /(lo*lo);
        // Vector3d s2 = -s1 - s3;
        // g(i*3*num_points_per_agent + 3*(j-1)+0) += K*s1(0); 
        // g(i*3*num_points_per_agent + 3*(j-1)+1) += K*s1(1); 
        // g(i*3*num_points_per_agent + 3*(j-1)+2) += K*s1(2); 
        // g(i*3*num_points_per_agent + 3*j+0)     += K*s2(0); 
        // g(i*3*num_points_per_agent + 3*j+1)     += K*s2(1); 
        // g(i*3*num_points_per_agent + 3*j+2)     += K*s2(2); 
        // g(i*3*num_points_per_agent + 3*(j+1)+0) += K*s3(0); 
        // g(i*3*num_points_per_agent + 3*(j+1)+1) += K*s3(1); 
        // g(i*3*num_points_per_agent + 3*(j+1)+2) += K*s3(2); 


        // double K = K_acc;
        // double t0 = -v1.dot(v2)*v1.dot(v2)/((v1.squaredNorm())*(v2.squaredNorm()));
        // long double t1 = -v1.dot(v2)*v1.dot(v2)/((v1.squaredNorm())*(v2.squaredNorm()))+1.0;
        // double t2 = ((q6-q9)*1.0/v1.norm()*1.0/v2.norm()-fabs(q3-q6)*(sign(q3,q6))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/pow(v1.squaredNorm(),3.0/2.0)*1.0/v2.norm());
        // long double t3 = sqrt(t1 + eps);
        // double t4 = (acos(v1.dot(v2)/(v1.norm()*v2.norm())));
        // double g1 = -K*1.0/t3*t2*t4;
        // std::cout<<"things"<<std::endl;
        // std::cout<<t0<<", "<<t1<<", "<<t2<<", "<<t3<<", "<<t4<<", "<<g1<<std::endl;
 

        // g(i*3*num_points_per_agent + 3*(j-1)+0) += -K*1.0/sqrt(-v1.dot(v2)*v1.dot(v2)/((v1.squaredNorm())*(v2.squaredNorm()))+1.0 + eps)*((q4-q7)*1.0/v1.norm()*1.0/v2.norm()-fabs(q1-q4)*(sign(q1,q4))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/pow(v1.squaredNorm(),3.0/2.0)*1.0/v2.norm())*(acos(v1.dot(v2)/(v1.norm()*v2.norm())));
        // g(i*3*num_points_per_agent + 3*(j-1)+1) += -K*1.0/sqrt(-v1.dot(v2)*v1.dot(v2)/((v1.squaredNorm())*(v2.squaredNorm()))+1.0 + eps)*((q5-q8)*1.0/v1.norm()*1.0/v2.norm()-fabs(q2-q5)*(sign(q2,q5))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/pow(v1.squaredNorm(),3.0/2.0)*1.0/v2.norm())*(acos(v1.dot(v2)/(v1.norm()*v2.norm())));
        // g(i*3*num_points_per_agent + 3*(j-1)+2) += -K*1.0/sqrt(-v1.dot(v2)*v1.dot(v2)/((v1.squaredNorm())*(v2.squaredNorm()))+1.0 + eps)*((q6-q9)*1.0/v1.norm()*1.0/v2.norm()-fabs(q3-q6)*(sign(q3,q6))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/pow(v1.squaredNorm(),3.0/2.0)*1.0/v2.norm())*(acos(v1.dot(v2)/(v1.norm()*v2.norm())));

        // g(i*3*num_points_per_agent + 3*j+0) += -K*1.0/sqrt(-v1.dot(v2)*v1.dot(v2)/((v1.squaredNorm())*(v2.squaredNorm()))+1.0 + eps)*(acos(v1.dot(v2)/(v1.norm()*v2.norm())))*((q1-q4*2.0+q7)*1.0/v1.norm()*1.0/v2.norm()+fabs(q1-q4)*(sign(q1,q4))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/pow(v1.squaredNorm(),3.0/2.0)*1.0/v2.norm()-fabs(q4-q7)*(sign(q4,q7))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/v1.norm()*1.0/pow(v2.squaredNorm(),3.0/2.0));
        // g(i*3*num_points_per_agent + 3*j+1) += -K*1.0/sqrt(-v1.dot(v2)*v1.dot(v2)/((v1.squaredNorm())*(v2.squaredNorm()))+1.0 + eps)*(acos(v1.dot(v2)/(v1.norm()*v2.norm())))*((q2-q5*2.0+q8)*1.0/v1.norm()*1.0/v2.norm()+fabs(q2-q5)*(sign(q2,q5))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/pow(v1.squaredNorm(),3.0/2.0)*1.0/v2.norm()-fabs(q5-q8)*(sign(q5,q8))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/v1.norm()*1.0/pow(v2.squaredNorm(),3.0/2.0));
        // g(i*3*num_points_per_agent + 3*j+2) += -K*1.0/sqrt(-v1.dot(v2)*v1.dot(v2)/((v1.squaredNorm())*(v2.squaredNorm()))+1.0 + eps)*(acos(v1.dot(v2)/(v1.norm()*v2.norm())))*((q3-q6*2.0+q9)*1.0/v1.norm()*1.0/v2.norm()+fabs(q3-q6)*(sign(q3,q6))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/pow(v1.squaredNorm(),3.0/2.0)*1.0/v2.norm()-fabs(q6-q9)*(sign(q6,q9))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/v1.norm()*1.0/pow(v2.squaredNorm(),3.0/2.0));
        
        // g(i*3*num_points_per_agent + 3*(j+1)+0) += K*1.0/sqrt(-v1.dot(v2)*v1.dot(v2)/((v1.squaredNorm())*(v2.squaredNorm()))+1.0 + eps)*((q1-q4)*1.0/v1.norm()*1.0/v2.norm()-fabs(q4-q7)*(sign(q4,q7))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/v1.norm()*1.0/pow(v2.squaredNorm(),3.0/2.0))*(acos(v1.dot(v2)/(v1.norm()*v2.norm())));
        // g(i*3*num_points_per_agent + 3*(j+1)+1) += K*1.0/sqrt(-v1.dot(v2)*v1.dot(v2)/((v1.squaredNorm())*(v2.squaredNorm()))+1.0 + eps)*((q2-q5)*1.0/v1.norm()*1.0/v2.norm()-fabs(q5-q8)*(sign(q5,q8))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/v1.norm()*1.0/pow(v2.squaredNorm(),3.0/2.0))*(acos(v1.dot(v2)/(v1.norm()*v2.norm())));
        // g(i*3*num_points_per_agent + 3*(j+1)+2) += K*1.0/sqrt(-v1.dot(v2)*v1.dot(v2)/((v1.squaredNorm())*(v2.squaredNorm()))+1.0 + eps)*((q3-q6)*1.0/v1.norm()*1.0/v2.norm()-fabs(q6-q9)*(sign(q6,q9))*((q1-q4)*(q4-q7)+(q2-q5)*(q5-q8)+(q3-q6)*(q6-q9))*1.0/v1.norm()*1.0/pow(v2.squaredNorm(),3.0/2.0))*(acos(v1.dot(v2)/(v1.norm()*v2.norm())));

        //----
        // if(isnan(g(i*3*num_points_per_agent + 3*(j-1)+0)))
        //   g(i*3*num_points_per_agent + 3*(j-1)+0) = 0;
        // if(isnan(g(i*3*num_points_per_agent + 3*(j-1)+1)))
        //   g(i*3*num_points_per_agent + 3*(j-1)+1) = 0;
        // if(isnan(g(i*3*num_points_per_agent + 3*(j-1)+2)))
        //   g(i*3*num_points_per_agent + 3*(j-1)+2) = 0;

        // if(isnan(g(i*3*num_points_per_agent + 3*(j)+0)))
        //   g(i*3*num_points_per_agent + 3*(j)+0) = 0;
        // if(isnan(g(i*3*num_points_per_agent + 3*(j)+1)))
        //   g(i*3*num_points_per_agent + 3*(j)+1) = 0;
        // if(isnan(g(i*3*num_points_per_agent + 3*(j)+2))){
        //   g(i*3*num_points_per_agent + 3*(j)+2) = 0;
        //   exit(0);
        // }

        // if(isnan(g(i*3*num_points_per_agent + 3*(j+1)+0)))
        //   g(i*3*num_points_per_agent + 3*(j+1)+0) = 0;
        // if(isnan(g(i*3*num_points_per_agent + 3*(j+1)+1)))
        //   g(i*3*num_points_per_agent + 3*(j+1)+1) = 0;
        // if(isnan(g(i*3*num_points_per_agent + 3*(j+1)+2)))
        //   g(i*3*num_points_per_agent + 3*(j+1)+2) = 0;
        //----

      }

    }
  }
}