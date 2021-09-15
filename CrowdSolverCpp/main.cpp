#include <igl/opengl/glfw/Viewer.h>
#include <cstdio>
#include <json.hpp>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include "kinetic_energy.h"
#include "kinetic_gradient.h"
#include "kinetic_hessian.h"
#include "tol_energy.h"
#include "tol_gradient.h"
#include "tol_hessian.h"
#include "accel_energy.h"
#include "accel_gradient.h"
#include "accel_hessian.h"

using namespace Eigen;
using namespace crowds;
using json = nlohmann::json;

VectorXd q,t, beq, b;
MatrixXd Aeq, A;//Todo: make sparse

int num_agents;
VectorXd UserTols;
//Weights
double K_agent = 1;
double K_tol =   1; 
double K_ke =    1;
double K_reg =   1;
double K_acc =   1;
std::vector<VectorXd> Qvec;



json j_input;
//std::string fname = "../../Scenes/output_results/complex_maze/square_maze/one_agent/";
std::string fname = "../../Scenes/output_results/scaling_tests/2_agents/";

void add_agent_equality_contraints(std::vector<int>& to_fix, VectorXd& av, MatrixXd& Aeq, VectorXd& beq)
{
  std::cout<<Aeq<<std::endl;
  MatrixXd A = MatrixXd::Zero(to_fix.size(), av.size()+1); //Extra +1 for the tol variable //Todo: make sparse,
  VectorXd b = VectorXd::Zero(to_fix.size());

  for(int i=0; i<to_fix.size(); i++)
  {
    b(i) = av(to_fix[i]);
    A(i, to_fix[i]) = 1;
  }

  MatrixXd newAeq = MatrixXd::Zero(Aeq.rows() + A.rows(), Aeq.cols() + A.cols());
  newAeq.block(0,0,Aeq.rows(), Aeq.cols()) = Aeq;
  newAeq.block(Aeq.rows(), Aeq.cols(), A.rows(), A.cols()) = A;

  VectorXd newbeq(beq.size()+b.size());
  newbeq<<beq,b;


  Aeq = newAeq;
  beq = newbeq;
}

void readInit()
{
  std::cout<<"reading file"<<std::endl;
  std::ifstream input_file(fname + "initial.json");
  input_file >> j_input;
  std::cout<<j_input["agents"][0].size()<<std::endl;
  std::cout<<"setting up agents"<<std::endl;
  //setup agent positions  
  for(int i=1; i<j_input["agents"].size(); ++i)
  {
    // VectorXd av = VectorXd::Zero(3*(j_input["agents"][i]["v"].size()));
    // //std::cout<<"av size: "<<j_input["agents"]["v"]<<std::endl;
    // for(int j=0; j<j_input["agents"][i]["v"].size(); j++)
    // {
    //   av.segment<3>(3*(j)) = Vector3d(j_input["agents"][i]["v"][j][0], 
    //                                 j_input["agents"][i]["v"][j][1], 
    //                                 j_input["agents"][i]["v"][j][2]);
    // }
    VectorXd av = VectorXd::Zero(3*3);
    av[0] = 0;
    av[1] = 0;
    av[2] = 0;
    av[3] = 0.5;
    av[4] = 0.3;
    av[5] = 0.5;
    av[6] = 1;
    av[7] = 1;
    av[8] = 1;
    // av[9] = 1.2;
    // av[10] = 1.7;
    // av[11] = 1.1;


    // av[0] = 0.5;
    // av[1] = 0.3;
    // av[2] = 0.5;
    // av[3] = 1;
    // av[4] = 1;
    // av[5] = 1;
    // av[6] = 1.2;
    // av[7] = 1.7;
    // av[8] = 1.1;


    Qvec.push_back(av);

    //constraints
    std::vector<int> to_fix;
    to_fix.push_back(0); //x0
    to_fix.push_back(1); //y0
    to_fix.push_back(2); //t0
    to_fix.push_back(av.size()-3); //xn
    to_fix.push_back(av.size()-2); //yn
    to_fix.push_back(av.size()-1); //tn Todo: make this inequality constraint

    add_agent_equality_contraints(to_fix, av, Aeq, beq);

  }
  q.resize(Qvec.size()*Qvec[0].size());
  for(int i=0;i<Qvec.size(); i++){
    q.segment(Qvec[0].size()*i, Qvec[i].size()) = Qvec[i];
  }
  std::cout<<"NUM Agents: "<<Qvec.size()<<std::endl;
  std::cout<<"NUM Points Per: "<<Qvec[0].size()<<std::endl;
  std::cout<<"NUM DOFS: "<<q.size()<<std::endl;
  return;
}

void writeResults()
{
  return;
}


double path_energy(VectorXd& x, VectorXd& UserTols, int num_agents, int num_points_per_agent, std::vector<VectorXd>& Q)
{
  VectorXd q = x.head(x.size() - num_agents);
  VectorXd t = x.tail(num_agents);

  double e = 0.0;
  double e_ke = kinetic_energy(q, num_agents, num_points_per_agent, K_ke);
  double e_tol = tol_energy(t, UserTols, K_tol);
  e = e_ke + e_tol;
  return e;
}

void path_gradient(VectorXd& x, VectorXd& UserTols, int num_agents, int num_points_per_agent, std::vector<VectorXd>& Q, VectorXd& g)
{
  VectorXd q = x.head(x.size() - num_agents);
  VectorXd t = x.tail(num_agents);

  
  g.resize(x.size()); 
  g.setZero();

  //individual gradients
  VectorXd g_ke = VectorXd::Zero(q.size());
  VectorXd g_tol = VectorXd::Zero(num_agents);

  kinetic_gradient(q, num_agents, num_points_per_agent, K_ke, g_ke);
  tol_gradient(t, UserTols, K_tol, g_tol);

  g.head(q.size()) += g_ke;
  g.tail(t.size()) += g_tol;
}

void path_hessian(VectorXd& x, VectorXd& UserTols, int num_agents, int num_points_per_agent, std::vector<VectorXd>& Q, MatrixXd& H)
{
  VectorXd q = x.head(x.size() - num_agents);
  VectorXd t = x.tail(num_agents);

  H.resize(x.size(), x.size());
  H.setZero();

  //individual gradients
  MatrixXd H_ke = MatrixXd::Zero(q.size(), q.size());
  MatrixXd H_tol = MatrixXd::Zero(num_agents, num_agents);

  kinetic_hessian(q, num_agents, num_points_per_agent, K_ke, H_ke);
  tol_hessian(t, UserTols, K_tol, H_tol);


  H.block(0,0,q.size(), q.size()) = H_ke;
  H.block(q.size() , q.size(), t.size(), t.size()) = H_tol;
}

void fd_check_gradient()
{
  //x, UserTols, num_agents, scene, e, surf_anim
  int num_agents = Qvec.size();

  int num_points_per_agent = Qvec[0].size()/3;
  VectorXd UserTols = 0.75*VectorXd::Ones(num_agents);

  double eps = 1e-5;
  double e0 = accel_energy(q, num_agents, num_points_per_agent, K_acc);
  std::cout<<"e0: "<<e0<<std::endl;
  VectorXd fdg = VectorXd::Zero(q.size());
  for(int i=0; i<fdg.size(); ++i)
  {
    VectorXd gl, gr;
    q(i) += eps;
    double er = accel_energy(q, num_agents, num_points_per_agent, K_acc);
    q(i) -= eps;

   

    double fd = (er - e0)/eps;
    fdg(i) = fd;
  }

  VectorXd g;
  accel_gradient(q, num_agents, num_points_per_agent, K_ke, g);

  std::cout<<g.transpose()<<std::endl;
  std::cout<<"---------------------"<<std::endl;
  std::cout<<fdg.transpose()<<std::endl;
  std::cout<<"#######################"<<std::endl;
  std::cout<<g.transpose() - fdg.transpose()<<std::endl;
}

void fd_check_hessian()
{
  //x, UserTols, num_agents, scene, e, surf_anim
  int num_agents = Qvec.size();
  int num_points_per_agent = Qvec[0].size()/3;
  VectorXd UserTols = 0.75*VectorXd::Ones(num_agents);

  double eps = 1e-4;
  VectorXd grad;
  accel_gradient(q, num_agents, num_points_per_agent, K_ke, grad);

  MatrixXd fdH = MatrixXd::Zero(q.size(), q.size());
  for(int i=0; i<grad.size(); ++i)
  {
    VectorXd gl, gr;
    q(i) += 0.5*eps;
    accel_gradient(q, num_agents, num_points_per_agent, K_ke, gr);
    q(i) -= eps;

    accel_gradient(q, num_agents, num_points_per_agent, K_ke, gl);
    q(i) += 0.5*eps;

    VectorXd fd = (gr - gl)/eps;
    fdH.row(i) = fd;
  }

  MatrixXd H;
  accel_hessian(q, num_agents, num_points_per_agent, K_ke, H);
  std::cout<<"--------------"<<std::endl;
  std::cout<<H<<std::endl;
  std::cout<<"--------------"<<std::endl;
  std::cout<<fdH<<std::endl;
  std::cout<<"-#########--"<<std::endl;
  std::cout<<fdH - H<<std::endl;
}

void line_search()
{

}

void newton_iteration()
{

}

void solve()
{

  //Constants
  int num_agents = Qvec.size();
  int num_points_per_agent = Qvec[0].size()/3;
  VectorXd UserTols = 0.75*VectorXd::Ones(num_agents);
  MatrixXd Z = MatrixXd::Zero(Aeq.rows(), Aeq.rows());

  //variables
  VectorXd x(q.size()+UserTols.size());
  x << q, 1e-8*VectorXd::Ones(num_agents);

  VectorXd grad;
  MatrixXd Hess;
  // double E = path_energy(x, UserTols, num_agents, num_points_per_agent, Qvec);
  // path_gradient(x, UserTols, num_agents, num_points_per_agent, Qvec, g);
  // path_hessian(x, UserTols, num_agents, num_points_per_agent, Qvec, H);

  int MAX_ITERS = 100;

  VectorXd x_n = x; //new x
  VectorXd lmb_n = VectorXd::Zero(Aeq.rows());// kkt lambda

  for(int iters = 0 ; iters < MAX_ITERS; iters++){
    //get initial energy
    double E_init = path_energy(x_n, UserTols, num_agents, num_points_per_agent, Qvec);
    //get force
    path_gradient(x_n, UserTols, num_agents, num_points_per_agent, Qvec, grad);
    //get stiffness
    path_hessian(x_n, UserTols, num_agents, num_points_per_agent, Qvec, Hess);


    //setup KKT solve: https://math.dartmouth.edu/~m126w18/pdf/part6.pdf
    MatrixXd KKT_A(Hess.rows() + Aeq.rows(), Hess.rows() + Aeq.rows());
    KKT_A<<Hess, Aeq.transpose(), 
          Aeq, Z;

    std::cout<<Hess.rows()<<std::endl;
    std::cout<<Aeq.rows()<<","<<Aeq.cols()<<std::endl;
    std::cout<<Z.rows()<<","<<Z.cols()<<std::endl;
    std::cout<<grad.size()<<std::endl;
    std::cout<<x_n.size()<<std::endl;
    std::cout<<beq.size()<<std::endl;


    VectorXd KKT_b(grad.size() + beq.size());
    VectorXd Ax_b = Aeq*x_n - beq;
    KKT_b<<grad,Ax_b;

    FullPivLU<MatrixXd> solver(KKT_A);
    Eigen::VectorXd kkt_dx = -solver.solve(KKT_b);

    VectorXd dx = kkt_dx.head(x_n.size());



    double step = 1.0;
    // //line search
    // double condition = 1e-7*grad.transpose()*dx;

    // double E0 = path_energy(x_n, UserTols, num_agents, num_points_per_agent, Qvec);

    // for(int lsits =0; lsits<50; lsits++){
    //   Eigen::VectorXd ls_x_n = x_n + step*dx;
    //   double newE = path_energy(ls_x_n, UserTols, num_agents, num_points_per_agent, Qvec);
      
    //   if(newE - E0 > step*condition){
    //     step  *= 0.5;
    //   }else{
    //     break;
    //   }
    // }
    // //end line search

    x_n += 1*dx;

    //get final energy
    double E_final = path_energy(x_n, UserTols, num_agents, num_points_per_agent, Qvec);

    if(grad.norm() < 1e-3)
    {
      break;
    }
    std::cout<<"STATS"<<std::endl;
    std::cout<< E_init<<std::endl;
    std::cout<< E_final<<std::endl;
    std::cout<< grad.norm()<<std::endl;
    std::cout<< step<<std::endl;
    std::cout<<std::endl;
  }
  x = x_n;
  q = x.head(q.size());

}


bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    viewer.data().clear();

    int num_agents = Qvec.size();
    int num_points_per_agent = Qvec[0].size()/3;
    int num_total_points = q.size()/3;

    for(int i=0; i<num_agents; i++)
    {
      double mass_i = 1.0;//hard coded for now

      VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
      MatrixXd Q_i = Map<MatrixXd>(q_i.data(), 3, q_i.size()/3).transpose();

      MatrixXd P1 = Q_i.block(0,0,Q_i.rows()-1, 3);
      MatrixXd P2 = Q_i.block(1,0,Q_i.rows()-1,3);
      viewer.data().add_edges(P1,P2,Eigen::RowVector3d(1,1,0));
    }


    // MatrixXd P1 = MatrixXd::Zero(num_total_points - num_agents, 3);
    // MatrixXd P2 = MatrixXd::Zero(num_total_points - num_agents, 3);

    // for(int i=0; i<num_agents; i++){
    //   for(int j=0; j<Qvec[i].size()/3 -1; j++){
    //     P1(num_points_per_agent*i + j, 0) = Q[i][3*j + 0];
    //     P1(num_points_per_agent*i + j, 1) = Q[i][3*j + 1];
    //     P1(num_points_per_agent*i + j, 2) = Q[i][3*j + 2];

    //     P2(num_points_per_agent*i + j, 0) = Q[i][3*(j+1) + 0];
    //     P2(num_points_per_agent*i + j, 1) = Q[i][3*(j+1) + 1];
    //     P2(num_points_per_agent*i + j, 2) = Q[i][3*(j+1) + 2];
    //   }
    // }
    
    // viewer.data().add_edges(P1,P2,Eigen::RowVector3d(1,1,0));
    return true;

  }

  if( key == '2')
  {

  }

  return false;
}

int main(int argc, char *argv[])
{

  readInit();
  // fd_check_gradient();
  fd_check_hessian();
  exit(0);
  solve();
  MatrixXd Q = Map<MatrixXd>(q.data(), 3, q.size()/3).transpose();
  std::cout<<Q<<std::endl;
  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  //viewer.data().set_mesh(V, F);
  //viewer.data().set_face_based(true);

  viewer.callback_key_down = &key_down;
  viewer.launch();

  return 1;
}
