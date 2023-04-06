//
// Created by Jemin Hwangbo on 2022/03/17.
//

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "exercise1_20233386.hpp"
#include "raisim/RaisimServer.hpp"


int main(int argc, char* argv[]) {
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();

  // aliengo
  auto aliengo = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/aliengo/aliengo.urdf");
  aliengo->setName("aliengo");
  server.focusOn(aliengo);

  // aliengo configuration
  Eigen::VectorXd jointNominalConfig(aliengo->getGeneralizedCoordinateDim());
  jointNominalConfig << 6, 6, 0.54, 1.0, 1.0, 0.0, 1.0, 0.1, 0.4, -0.95, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8;
  aliengo->setGeneralizedCoordinate(jointNominalConfig);

  // debug sphere
  auto debugSphere = server.addVisualSphere("debug_sphere", 0.02);
  debugSphere->setColor(1,0,0,1);
  debugSphere->setPosition(getEndEffectorPosition(jointNominalConfig));

  // solution sphere
  auto answerSphere = server.addVisualSphere("answer_sphere", 0.02);
  answerSphere->setColor(0,1,0,1);
  raisim::Vec<3> pos;
  aliengo->getFramePosition("FR_foot_fixed", pos);
  answerSphere->setPosition(pos.e());

  if (pos.e() != getEndEffectorPosition(jointNominalConfig))
      std::cout << pos.e() << std::endl;
      std::cout <<getEndEffectorPosition(jointNominalConfig)<< std::endl;
  // visualization
  server.launchServer();
  for (int i=0; i<2000000; i++)
    std::this_thread::sleep_for(std::chrono::microseconds(1000));

  server.killServer();
}
