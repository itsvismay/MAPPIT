#include <igl/opengl/glfw/Viewer.h>
#include <cstdio>
#include <json.hpp>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include "kinetic_energy.h"
#include "kinetic_gradient.h"
#include "kinetic_hessian.h"
#include "reg_energy.h"
#include "reg_gradient.h"
#include "reg_hessian.h"
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
std::string fname = "../../Scenes/2_output_results/scaling_tests/2_agents/run12/";

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
  for(int i=0; i<j_input["agents"].size(); ++i)
  {
    VectorXd av = VectorXd::Zero(3*(j_input["agents"][i]["v"].size()));
    //std::cout<<"av size: "<<j_input["agents"]["v"]<<std::endl;
    for(int j=0; j<j_input["agents"][i]["v"].size(); j++)
    {
      av.segment<3>(3*(j)) = Vector3d(j_input["agents"][i]["v"][j][0], 
                                    j_input["agents"][i]["v"][j][1], 
                                    j_input["agents"][i]["v"][j][2]);
    }
  
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

  // for(int i=0; i<1; ++i)
  // {


  //   VectorXd av = VectorXd::Zero(3*3);
  //   av[0] = 0;
  //   av[1] = 0;
  //   av[2] = 0;
  //   av[3] = 0.5;
  //   av[4] = 0.3;
  //   av[5] = 0.5;
  //   av[6] = 1;
  //   av[7] = 1;
  //   av[8] = 1;
  //   // av[9] = 1.2;
  //   // av[10] = 1.7;
  //   // av[11] = 1.1;


  //   // av[0] =        28.125;
  //   // av[1] =            0;
  //   // av[2] =        0.005;
  //   // av[3] =           25;
  //   // av[4] =            0;
  //   // av[5] =        0.006;
  //   // av[6] =       21.875;
  //   // av[7] =            0;
  //   // av[8] =        0.007;
  
  //   Qvec.push_back(av);

  //   //constraints
  //   std::vector<int> to_fix;
  //   to_fix.push_back(0); //x0
  //   to_fix.push_back(1); //y0
  //   to_fix.push_back(2); //t0
  //   to_fix.push_back(av.size()-3); //xn
  //   to_fix.push_back(av.size()-2); //yn
  //   to_fix.push_back(av.size()-1); //tn Todo: make this inequality constraint

  //   add_agent_equality_contraints(to_fix, av, Aeq, beq);
  // }

  q.resize(Qvec.size()*Qvec[0].size());
  for(int i=0;i<Qvec.size(); i++){
    q.segment(Qvec[0].size()*i, Qvec[i].size()) = Qvec[i];
  }
  std::cout<<"NUM Agents: "<<Qvec.size()<<std::endl;
  std::cout<<"NUM Points Per: "<<Qvec[0].size()<<std::endl;
  std::cout<<"NUM DOFS: "<<q.size()<<std::endl;
  return;
}

void fd_check_gradient()
{
  //x, UserTols, num_agents, scene, e, surf_anim
  int num_agents = Qvec.size();

  int num_points_per_agent = Qvec[0].size()/3;

  double eps = 1e-5;
  double e0 = reg_energy(q, num_agents, num_points_per_agent, K_acc);
  std::cout<<"e0: "<<e0<<std::endl;
  VectorXd fdg = VectorXd::Zero(q.size());
  for(int i=0; i<fdg.size(); ++i)
  {
    q(i) += eps;
    double er = reg_energy(q, num_agents, num_points_per_agent, K_acc);
    q(i) -= eps;
    double fd = (er - e0)/eps;
    fdg(i) = fd;
  }

  VectorXd g;
  reg_gradient(q, num_agents, num_points_per_agent, K_ke, g);

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
  reg_gradient(q, num_agents, num_points_per_agent, K_ke, grad);

  MatrixXd fdH = MatrixXd::Zero(q.size(), q.size());
  for(int i=0; i<grad.size(); ++i)
  {
    VectorXd gl, gr;
    q(i) += 0.5*eps;
    reg_gradient(q, num_agents, num_points_per_agent, K_ke, gr);
    q(i) -= eps;

    reg_gradient(q, num_agents, num_points_per_agent, K_ke, gl);
    q(i) += 0.5*eps;

    VectorXd fd = (gr - gl)/eps;
    fdH.row(i) = fd;
  }

  SparseMatrix<double> H;
  reg_hessian(q, num_agents, num_points_per_agent, K_ke, H);
  MatrixXd fullH = MatrixXd(H);
  std::cout<<"--------------"<<std::endl;
  std::cout<<fullH.block<9,9>(0,0)<<std::endl;
  std::cout<<"--------------"<<std::endl;
  std::cout<<fdH.block<9,9>(0,0)<<std::endl;
  std::cout<<"-#########--"<<std::endl;
  MatrixXd Diff = fdH - fullH;
  std::cout<<Diff.block<9,9>(0,0)<<std::endl;
}


int main(int argc, char *argv[])
{

  readInit();
  //fd_check_gradient();
  fd_check_hessian();
  exit(0);

  return 1;
}
