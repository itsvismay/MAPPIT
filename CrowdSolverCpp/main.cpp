#include <igl/opengl/glfw/Viewer.h>
#include <cstdio>
#include <json.hpp>

using namespace Eigen;
using json = nlohmann::json;

VectorXd q;
int num_agents;
VectorXd UserTols;
//Weights
double K_agent = 1;
double K_tol =   1; 
double K_ke =    1;
double K_reg =   1;
std::vector<VectorXd> Qvec;



json j_input;
std::string fname = "../../Scenes/output_results/scaling_tests/1_agents/";

void readInit()
{
  std::ifstream input_file(fname + "initial.json");
  input_file >> j_input;
  std::cout<<j_input["agents"][0].size()<<std::endl;

  
  for(int i=0; i<j_input["agents"].size(); ++i)
  {
    VectorXd av = VectorXd::Zero(3*j_input["agents"][i]["v"].size());
    //std::cout<<"av size: "<<j_input["agents"]["v"]<<std::endl;
    for(int j=0; j<j_input["agents"][i]["v"].size(); j++)
    {
      av.segment<3>(3*j) = Vector3d(j_input["agents"][i]["v"][j][0], 
                                    j_input["agents"][i]["v"][j][1], 
                                    j_input["agents"][i]["v"][j][2]);
    }
    Qvec.push_back(av);
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

double kinetic_energy(VectorXd& q, int num_agents, int num_points_per_agent, double K_Ke)
{
  double e = 0.0;
  for(int i=0; i<num_agents; i++)
  {
    double mass_i = 1.0;//hard coded for now

    VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
    MatrixXd Q_i = Map<MatrixXd>(q_i.data(), 3, q_i.size()/3).transpose();

    VectorXd dx = Q_i.col(0).tail(num_points_per_agent -1) - Q_i.col(0).head(num_points_per_agent - 1) ;
    VectorXd dy = Q_i.col(1).tail(num_points_per_agent -1) - Q_i.col(1).head(num_points_per_agent - 1) ;
    VectorXd dt = Q_i.col(2).tail(num_points_per_agent -1) - Q_i.col(2).head(num_points_per_agent - 1) ;
    

    ArrayXd dx2dy2_dt = (dx.array()*dx.array() + dy.array()*dy.array())/dt.array();

    double e_i = 0.5*mass_i*dx2dy2_dt.matrix().sum();//e = 0.5*(dx^2 + dy^2)/dt
    e += e_i;
  }
}
void kinetic_gradient(VectorXd& q, int num_agents, int num_points_per_agent, double K_Ke, VectorXd& g)
{
  g.resize(q.size());
  g.setZero();

  double e = 0.0;
  for(int i=0; i<num_agents; i++)
  {
    double mass_i = 1.0;//hard coded for now

    VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
    
    //E = ((x2 - x1)^2 + (y2 - y1)^2)/(t2 - t1) -> dE/dx

    //For each rod segment, do this:
    for(int e =0; e<num_points_per_agent - 1; e++)
    {
      double x1 = q_i(3*e+0);
      double y1 = q_i(3*e+1);
      double t1 = q_i(3*e+2);
      double x2 = q_i(3*(e+1)+0);
      double y2 = q_i(3*(e+1)+1);
      double t2 = q_i(3*(e+1)+2);

      g(i*3*num_points_per_agent + 3*e+0) += -mass_i*(x2 - x1)/(t2 - t1); //dE/dx1
      g(i*3*num_points_per_agent + 3*e+1) += -mass_i*(y2 - y1)/(t2 - t1); //dE/dy1
      g(i*3*num_points_per_agent + 3*e+2) += 0.5*mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) )/((t2 - t1)*(t2 - t1)); //dE/dt1

      g(i*3*num_points_per_agent + 3*(e+1)+0) += mass_i*(x2 - x1)/(t2 - t1); //dE/dx2
      g(i*3*num_points_per_agent + 3*(e+1)+1) += mass_i*(y2 - y1)/(t2 - t1); //dE/dy2
      g(i*3*num_points_per_agent + 3*(e+1)+2) += -0.5*mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) )/((t2 - t1)*(t2 - t1)); //dE/dt2
    }

  }
}

void kinetic_hessian(VectorXd& q, int num_agents, int num_points_per_agent, double K_Ke, MatrixXd& H)
{
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
      H(i*3*num_points_per_agent + 3*e+0,     i*3*num_points_per_agent + 3*e+0)     += mass_i/(t2 - t1); 
      H(i*3*num_points_per_agent + 3*(e+1)+0, i*3*num_points_per_agent + 3*(e+1)+0) += mass_i/(t2 - t1);

      //ddE/dy1dy1, //ddE/dy2dy2
      H(i*3*num_points_per_agent + 3*e+1,     i*3*num_points_per_agent + 3*e+1)     += mass_i/(t2 - t1); 
      H(i*3*num_points_per_agent + 3*(e+1)+1, i*3*num_points_per_agent + 3*(e+1)+1) += mass_i/(t2 - t1); 

      //ddE/dx1dx2, ddE/dy1dy2
      H(i*3*num_points_per_agent + 3*e+0,     i*3*num_points_per_agent + 3*(e+1)+0)     += -mass_i/(t2 - t1); 
      H(i*3*num_points_per_agent + 3*(e+1)+0,     i*3*num_points_per_agent + 3*e+0)     += -mass_i/(t2 - t1); 

      H(i*3*num_points_per_agent + 3*e+1,     i*3*num_points_per_agent + 3*(e+1)+1)     += -mass_i/(t2 - t1); 
      H(i*3*num_points_per_agent + 3*(e+1)+1,     i*3*num_points_per_agent + 3*e+1)     += -mass_i/(t2 - t1); 

      //-----------------
      //ddE/dt1dx1, ddE/dx1dt1
      H(i*3*num_points_per_agent + 3*e+0,     i*3*num_points_per_agent + 3*e+2)     += -mass_i*(x2 - x1)/((t2-t1)*(t2 - t1));
      H(i*3*num_points_per_agent + 3*e+2,     i*3*num_points_per_agent + 3*e+0)     += -mass_i*(x2 - x1)/((t2-t1)*(t2 - t1));

      //ddE/dt1dx2, ddE/dx2dt1
      H(i*3*num_points_per_agent + 3*(e+1)+0, i*3*num_points_per_agent + 3*e+2)     += mass_i*(x2 - x1)/((t2-t1)*(t2 - t1));
      H(i*3*num_points_per_agent + 3*e+2,     i*3*num_points_per_agent + 3*(e+1)+0) += mass_i*(x2 - x1)/((t2-t1)*(t2 - t1));

      //ddE/dt2dx1, ddE/dx1dt2
      H(i*3*num_points_per_agent + 3*e+0,     i*3*num_points_per_agent + 3*(e+1)+2) += mass_i*(x2 - x1)/((t2-t1)*(t2 - t1));
      H(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*e+0)     += mass_i*(x2 - x1)/((t2-t1)*(t2 - t1));

      //ddE/dt2dx2, ddE/dx2dt2
      H(i*3*num_points_per_agent + 3*(e+1)+0, i*3*num_points_per_agent + 3*(e+1)+2) += -mass_i*(x2 - x1)/((t2-t1)*(t2 - t1));
      H(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*(e+1)+0) += -mass_i*(x2 - x1)/((t2-t1)*(t2 - t1));

      //-----------------
      //ddE/dt1dy1, ddE/dy1dt1
      H(i*3*num_points_per_agent + 3*e+1,     i*3*num_points_per_agent + 3*e+2)     += -mass_i*(y2 - y1)/((t2-t1)*(t2 - t1));
      H(i*3*num_points_per_agent + 3*e+2,     i*3*num_points_per_agent + 3*e+1)     += -mass_i*(y2 - y1)/((t2-t1)*(t2 - t1));

      //ddE/dt1dy2, ddE/dy2dt1
      H(i*3*num_points_per_agent + 3*(e+1)+1, i*3*num_points_per_agent + 3*e+2)     += mass_i*(y2 - y1)/((t2-t1)*(t2 - t1));
      H(i*3*num_points_per_agent + 3*e+2,     i*3*num_points_per_agent + 3*(e+1)+1) += mass_i*(y2 - y1)/((t2-t1)*(t2 - t1));

      //ddE/dt2dy1, ddE/dy1dt2
      H(i*3*num_points_per_agent + 3*e+1,     i*3*num_points_per_agent + 3*(e+1)+2) += -mass_i*(y2 - y1)/((t2-t1)*(t2 - t1));
      H(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*e+1)     += -mass_i*(y2 - y1)/((t2-t1)*(t2 - t1));

      //ddE/dt2dy2, ddE/dy2dt2
      H(i*3*num_points_per_agent + 3*(e+1)+1, i*3*num_points_per_agent + 3*(e+1)+2) += mass_i*(y2 - y1)/((t2-t1)*(t2 - t1));
      H(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*(e+1)+1) += mass_i*(y2 - y1)/((t2-t1)*(t2 - t1));

      //-----------------
      //ddE/dt1dt1
      H(i*3*num_points_per_agent + 3*e+2,     i*3*num_points_per_agent + 3*e+2)     += mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1))/((t2 - t1)*(t2 - t1)*(t2 - t1));

      //ddE/dt1dt2, ddE/dt1dt2 
      H(i*3*num_points_per_agent + 3*e+2,     i*3*num_points_per_agent + 3*(e+1)+2) += -mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1))/((t2 - t1)*(t2 - t1)*(t2 - t1));
      H(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*e+2)     += -mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1))/((t2 - t1)*(t2 - t1)*(t2 - t1));

      //ddE/dt2dt2
      H(i*3*num_points_per_agent + 3*(e+1)+2, i*3*num_points_per_agent + 3*(e+1)+2) += mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1))/((t2 - t1)*(t2 - t1)*(t2 - t1));

      //GRADIENTS
      // g(i*3*num_points_per_agent + 3*e+0) += -mass_i*(x2 - x1)/(t2 - t1); //dE/dx1
      // g(i*3*num_points_per_agent + 3*e+1) += -mass_i*(y2 - y1)/(t2 - t1); //dE/dy1
      // g(i*3*num_points_per_agent + 3*e+2) += 0.5*mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) )/((t2 - t1)*(t2 - t1)); //dE/dt1

      // g(i*3*num_points_per_agent + 3*(e+1)+0) += mass_i*(x2 - x1)/(t2 - t1); //dE/dx2
      // g(i*3*num_points_per_agent + 3*(e+1)+1) += mass_i*(y2 - y1)/(t2 - t1); //dE/dy2
      // g(i*3*num_points_per_agent + 3*(e+1)+2) += -0.5*mass_i*((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) )/((t2 - t1)*(t2 - t1)); //dE/dt2
    }

  }
}

double path_energy(VectorXd& x, VectorXd& UserTols, int num_agents, int num_points_per_agent, std::vector<VectorXd>& Q)
{
  VectorXd q = x.head(x.size() - num_agents);
  VectorXd T = x.tail(num_agents);


  double e_ke = kinetic_energy(q, num_agents, num_points_per_agent, K_ke);

}

void path_gradient(VectorXd& x, VectorXd& UserTols, int num_agents, int num_points_per_agent, std::vector<VectorXd>& Q, VectorXd& g)
{
  VectorXd q = x.head(x.size() - num_agents);
  VectorXd T = x.tail(num_agents);

  
  g.resize(x.size()); 
  g.setZero();

  //individual gradients
  VectorXd g_ke =VectorXd::Zero(q.size());

  kinetic_gradient(q, num_agents, num_points_per_agent, K_ke, g_ke);
  g.head(q.size()) += g_ke;

}

void path_hessian(VectorXd& x, VectorXd& UserTols, int num_agents, int num_points_per_agent, std::vector<VectorXd>& Q, MatrixXd& H)
{
  VectorXd q = x.head(x.size() - num_agents);
  VectorXd T = x.tail(num_agents);

  H.resize(x.size(), x.size());
  H.setZero();

  //individual gradients
  MatrixXd H_ke = MatrixXd::Zero(q.size(), q.size());

  kinetic_hessian(q, num_agents, num_points_per_agent, K_ke, H_ke);
  std::cout<<H_ke<<std::endl;
  
}

void fd_check()
{
  //x, UserTols, num_agents, scene, e, surf_anim
  int num_agents = Qvec.size();
  int num_points_per_agent = Qvec[0].size()/3;
  VectorXd UserTols = 0.75*VectorXd::Ones(num_agents);

  // double E = path_energy(x, UserTols, num_agents, num_points_per_agent, Qvec);

  // VectorXd g;
  // path_gradient(x, UserTols, num_agents, num_points_per_agent, Qvec, g);

  // MatrixXd H;
  // path_hessian(x, UserTols, num_agents, num_points_per_agent, Qvec, H);

  double eps = 1e-5;
  VectorXd grad;
  kinetic_gradient(q, num_agents, num_points_per_agent, K_ke, grad);

  MatrixXd fdH = MatrixXd::Zero(q.size(), q.size());
  for(int i=0; i<grad.size(); ++i)
  {
    VectorXd gl, gr;
    q(i) += 0.5*eps;
    kinetic_gradient(q, num_agents, num_points_per_agent, K_ke, gr);
    q(i) -= eps;

    kinetic_gradient(q, num_agents, num_points_per_agent, K_ke, gl);
    q(i) += 0.5*eps;

    VectorXd fd = (gr - gl)/eps;
    fdH.row(i) = fd;
  }

  MatrixXd H;
  kinetic_hessian(q, num_agents, num_points_per_agent, K_ke, H);
 

}

void solve()
{

  fd_check();
  exit(0);
  //x, UserTols, num_agents, scene, e, surf_anim
  int num_agents = Qvec.size();
  int num_points_per_agent = Qvec[0].size()/3;
  VectorXd UserTols = 0.75*VectorXd::Ones(num_agents);
  
  VectorXd x(q.size()+UserTols.size());
  x << q, 1e-8*VectorXd::Ones(num_agents);

  double E = path_energy(x, UserTols, num_agents, num_points_per_agent, Qvec);

  VectorXd g;
  path_gradient(x, UserTols, num_agents, num_points_per_agent, Qvec, g);

  MatrixXd H;
  path_hessian(x, UserTols, num_agents, num_points_per_agent, Qvec, H);


}


bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    viewer.data().clear();

    int num_agents = Qvec.size();
    int num_points_per_agent = Qvec[0].size()/3;
    int num_total_points = q.size()/3;

    MatrixXd P1 = MatrixXd::Zero(num_total_points - num_agents, 3);
    MatrixXd P2 = MatrixXd::Zero(num_total_points - num_agents, 3);

    for(int i=0; i<num_agents; i++){
      for(int j=0; j<Qvec[i].size()/3 -1; j++){
        std::cout<<num_points_per_agent*i + j<<std::endl;
        P1(num_points_per_agent*i + j, 0) = Qvec[i][3*j + 0];
        P1(num_points_per_agent*i + j, 1) = Qvec[i][3*j + 1];
        P1(num_points_per_agent*i + j, 2) = Qvec[i][3*j + 2];

        P2(num_points_per_agent*i + j, 0) = Qvec[i][3*(j+1) + 0];
        P2(num_points_per_agent*i + j, 1) = Qvec[i][3*(j+1) + 1];
        P2(num_points_per_agent*i + j, 2) = Qvec[i][3*(j+1) + 2];
      }
    }
    
    viewer.data().add_edges(P1,P2,Eigen::RowVector3d(1,1,0));
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
  solve();
  exit(0);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  //viewer.data().set_mesh(V, F);
  //viewer.data().set_face_based(true);

  viewer.callback_key_down = &key_down;
  viewer.launch();

  return 1;
}
