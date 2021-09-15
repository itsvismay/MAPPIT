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


      for(int j=1; j<Q_i.rows() -1; j++){
        //for each node thats not the first node (0) and last node
        double q1 = Q_i.row(j-1)[0];
        double q2 = Q_i.row(j-1)[1];
        double q3 = Q_i.row(j-1)[2];
        double q4 = Q_i.row(j)[0];
        double q5 = Q_i.row(j)[1];
        double q6 = Q_i.row(j)[2];
        double q7 = Q_i.row(j+1)[0];
        double q8 = Q_i.row(j+1)[1];
        double q9 = Q_i.row(j+1)[2];

        Vector3d v1 = Q_i.row(j) - Q_i.row(j-1);
        Vector3d v2 = Q_i.row(j+1) - Q_i.row(j);
        Vector3d v1xv2 = (v1.cross(v2)).transpose();
        double v1dv2 = v1.dot(v2);
        double v1xv2norm = v1xv2.norm();

        // //----Dave's Energy-----
        double eps = 1e-5;
        if(v1dv2<=eps){
          v1dv2 = eps;
          v2(0) += eps;
          v2(1) += eps;
          v2(2) += eps;
        }
       
        v1xv2 = (v1.cross(v2)).transpose();
        v1dv2 = v1.dot(v2);
        v1xv2norm = v1xv2.norm();
        // std::cout<<"v1: "<<v1.transpose()<<std::endl;
        // std::cout<<"v2: "<<v2.transpose()<<std::endl;
        // std::cout<<"v1dv2: "<<v1dv2<<std::endl;
        // std::cout<<"v1xv2: "<<v1xv2.transpose()<<std::endl;
        // std::cout<<"v1dv2norm: "<<v1xv2norm<<std::endl;

        double X = v1dv2;
        double Y = v1xv2norm;

        double dE = 0;
        double datan = 1.0/(1.0 + (Y/X)*(Y/X));
        if(v1xv2.norm()<=eps)
        {
          dE = 0;
        }
        else
        {
          dE = K_acc*atan2(Y, X);
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
        
        // std::cout<<"lodhi1: "<<lodhi1.transpose()<<std::endl;
        // std::cout<<"hidlo1: "<<hidlo1.transpose()<<std::endl;
        // std::cout<<"lodhi1 - hidlo1: "<<lodhi1.transpose() - hidlo1.transpose()<<std::endl;
        // std::cout<<"lolo: "<<lolo<<std::endl;
        // std::cout<<"dE: "<<dE<<std::endl;
   
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
        
        //------------------Symbolic gives me nans

        // double v11 = v1(0);
        // double v12 = v1(1);
        // double v13 = v1(2);

        // double v21 = v2(0);
        // double v22 = v2(1);
        // double v23 = v2(2);

        // Vector3d T1, T2;
        // T1(0) = -v21*1.0/pow(v1dv2,2.0)*sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0))+((v22*fabs(v11*v22-v12*v21)*(((v11*v22-v12*v21)/fabs(v11*v22-v12*v21)))*2.0+v23*fabs(v11*v23-v13*v21)*(((v11*v23-v13*v21)/fabs(v11*v23-v13*v21)))*2.0)*1.0/sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0)))/(v11*v21*2.0+v12*v22*2.0+v13*v23*2.0);
        // T1(1) = -v22*1.0/pow(v1dv2,2.0)*sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0))-((v21*fabs(v11*v22-v12*v21)*(((v11*v22-v12*v21)/fabs(v11*v22-v12*v21)))*2.0-v23*fabs(v12*v23-v13*v22)*(((v12*v23-v13*v22)/fabs(v12*v23-v13*v22)))*2.0)*1.0/sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0)))/(v11*v21*2.0+v12*v22*2.0+v13*v23*2.0);
        // T1(2) = -v23*1.0/pow(v1dv2,2.0)*sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0))-((v21*fabs(v11*v23-v13*v21)*(((v11*v23-v13*v21)/fabs(v11*v23-v13*v21)))*2.0+v22*fabs(v12*v23-v13*v22)*(((v12*v23-v13*v22)/fabs(v12*v23-v13*v22)))*2.0)*1.0/sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0)))/(v11*v21*2.0+v12*v22*2.0+v13*v23*2.0);

        // T2(0) = -v11*1.0/pow(v11*v21+v12*v22+v13*v23,2.0)*sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0))-((v12*fabs(v11*v22-v12*v21)*(((v11*v22-v12*v21)/fabs(v11*v22-v12*v21)))*2.0+v13*fabs(v11*v23-v13*v21)*(((v11*v23-v13*v21)/fabs(v11*v23-v13*v21)))*2.0)*1.0/sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0)))/(v11*v21*2.0+v12*v22*2.0+v13*v23*2.0);
        // T2(1) = -v12*1.0/pow(v11*v21+v12*v22+v13*v23,2.0)*sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0))+((v11*fabs(v11*v22-v12*v21)*(((v11*v22-v12*v21)/fabs(v11*v22-v12*v21)))*2.0-v13*fabs(v12*v23-v13*v22)*(((v12*v23-v13*v22)/fabs(v12*v23-v13*v22)))*2.0)*1.0/sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0)))/(v11*v21*2.0+v12*v22*2.0+v13*v23*2.0);
        // T2(2) = -v13*1.0/pow(v11*v21+v12*v22+v13*v23,2.0)*sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0))+((v11*fabs(v11*v23-v13*v21)*(((v11*v23-v13*v21)/fabs(v11*v23-v13*v21)))*2.0+v12*fabs(v12*v23-v13*v22)*(((v12*v23-v13*v22)/fabs(v12*v23-v13*v22)))*2.0)*1.0/sqrt(pow(fabs(v11*v22-v12*v21),2.0)+pow(fabs(v11*v23-v13*v21),2.0)+pow(fabs(v12*v23-v13*v22),2.0)))/(v11*v21*2.0+v12*v22*2.0+v13*v23*2.0);

        // std::cout<<lhs*-T1<<std::endl;
        // std::cout<<lhs*T2<<std::endl;
        // //Vector3d Gjprev = lhs * -T1;
        // //Vector3d Gjnext = lhs * T2;
        // //Vector3d Gj =  - Gjprev - Gjnext;
        //------------------

        g(i*3*num_points_per_agent + 3*(j-1)+0) += Gjprev(0);
        g(i*3*num_points_per_agent + 3*(j-1)+1) += Gjprev(1);
        g(i*3*num_points_per_agent + 3*(j-1)+2) += Gjprev(2);

        g(i*3*num_points_per_agent + 3*j+0)     += Gj(0);
        g(i*3*num_points_per_agent + 3*j+1)     += Gj(1);
        g(i*3*num_points_per_agent + 3*j+2)     += Gj(2);

        g(i*3*num_points_per_agent + 3*(j+1)+0) += Gjnext(0);
        g(i*3*num_points_per_agent + 3*(j+1)+1) += Gjnext(1);
        g(i*3*num_points_per_agent + 3*(j+1)+2) += Gjnext(2);

        //-------Etienne's Energy-----------------------------------
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

        // g(i*3*num_points_per_agent + 3*(j-1)+0) += Gjprev(0);
        // g(i*3*num_points_per_agent + 3*(j-1)+1) += Gjprev(1);
        // g(i*3*num_points_per_agent + 3*(j-1)+2) += Gjprev(2);

        // g(i*3*num_points_per_agent + 3*j+0)     += Gj(0);
        // g(i*3*num_points_per_agent + 3*j+1)     += Gj(1);
        // g(i*3*num_points_per_agent + 3*j+2)     += Gj(2);

        // g(i*3*num_points_per_agent + 3*(j+1)+0) += Gjnext(0);
        // g(i*3*num_points_per_agent + 3*(j+1)+1) += Gjnext(1);
        // g(i*3*num_points_per_agent + 3*(j+1)+2) += Gjnext(2);


      }

    }
  }
}