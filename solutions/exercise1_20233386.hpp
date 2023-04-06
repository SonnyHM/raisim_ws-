//
// Created by Jemin Hwangbo on 2022/03/17.
//

#ifndef ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_
#define ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_

#include <Eigen/Core>
#include "Eigen/Geometry"

/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {
  //////////////////////////
  ///// Your Code Here /////
  /////////////////////////
  Eigen::Quaterniond b_q(gc(3),gc(4),gc(5),gc(6));
  Eigen::Vector3d end_pos;
  Eigen::Vector3d r_0, r_1, r_2, r_3, r_4;
  Eigen::Matrix3d Rot_0, Rot_1, Rot_2, Rot_3;
  double th_1 = gc(7), th_2 = gc(8), th_3 = gc(9);

  //transpose quarternion of base to rotation matrix
  Rot_0 = b_q.normalized().toRotationMatrix();

  Rot_1 << 1, 0, 0,
           0, cos(th_1), -sin(th_1),
           0, sin(th_1), cos(th_1);

  Rot_2 << cos(th_2), 0, sin(th_2),
           0, 1, 0,
           -sin(th_2), 0, cos(th_2);

  Rot_3 << cos(th_3), 0, sin(th_3),
           0, 1, 0,
           -sin(th_3), 0, cos(th_3);

  r_4 << 0, 0, -0.25;
  r_3 << 0, 0, -0.25;
  r_2 << 0, -0.083, 0;
  r_1 << 0.2399, -0.051, 0;
  r_0 << gc(0), gc(1), gc(2);

//other rotation matrices are all I
  end_pos = r_0 + Rot_0*(r_1 + Rot_1*(r_2 + Rot_2*(r_3 + Rot_3 * r_4)));

  return end_pos; /// replace this
}

#endif // ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_
