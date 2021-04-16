#include <igl/opengl/glfw/Viewer.h>
#include <cstdio>
#include <json.hpp>

using namespace Eigen;
using json = nlohmann::json;

struct Particle
{
  int id;
  Vector3d pos;
  Vector3d prevPos;
  Vector3d vel;
  bool fixed = false;
  int nextId = -1;//id of the next particle in the chain
  int prevId = -1;//id of the previous particle in the chain
};

json j_input;
double dt = 1.0/60.0;
double radius = 0.25;
double substeps = 20;
std::vector<Particle> particles;
std::string fname = "../../Scenes/output_results/scaling_tests/test/";

void readInit()
{
  std::ifstream input_file(fname + "initial.json");
  input_file >> j_input;
  //std::cout<<j_input["agents"]<<std::endl;

  for(int i=0; i<j_input["agents"]["v"].size(); ++i){
    Particle p;
    p.id = i;
    p.pos = Vector3d(j_input["agents"]["v"][i][0], j_input["agents"]["v"][i][1], j_input["agents"]["v"][i][2]);
    p.prevPos = p.pos;
    p.vel.setZero();

    if(i==0)
    {  
      p.fixed = true;
    }
    if(i!=j_input["agents"]["v"].size()-1)
    {  
      p.nextId = i+1;
      p.fixed = true;
    }

    particles.push_back(p);
  }
  return;
}

void writeResults()
{
  return;
}

void simulate()
{

  double h = dt/substeps;

  for(int steps=0; steps<(int)substeps; steps ++)
  {
    //predict
    for(int i=0; i< particles.size(); i++)
    {
      particles[i].prevPos = particles[i].pos;
      particles[i].pos = h*particles[i].vel;
    }

    //update collisions

    //
  }

  
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    viewer.data().clear();

    MatrixXd P1 = MatrixXd::Zero(particles.size()-1, 3);
    MatrixXd P2 = MatrixXd::Zero(particles.size()-1, 3);
    for(int i=0; i<particles.size()-1; i++)
    {
      P1.row(i) = particles[i].pos;
      P2.row(i) = particles[i+1].pos;
    }

    P1 = P1/100.0;
    P2 = P2/100.0;
    viewer.data().add_edges(P1,P2,Eigen::RowVector3d(1,1,0));

    return true;

  }
  return false;
}

int main(int argc, char *argv[])
{

  readInit();

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  //viewer.data().set_mesh(V, F);
  //viewer.data().set_face_based(true);

  viewer.callback_key_down = &key_down;
  viewer.launch();

  return 1;
}
