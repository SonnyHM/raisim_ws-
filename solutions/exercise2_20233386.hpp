#pragma once
#include <Eigen/Core>
#include "Eigen/Geometry"

inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc, int i) {

    Eigen::Quaterniond b_q(gc(3),gc(4),gc(5),gc(6));
    Eigen::Vector3d end_pos;
    Eigen::Vector3d r_0, r_1, r_2, r_3, r_4;
    Eigen::Matrix3d Rot_0, Rot_1, Rot_2, Rot_3;
    double th_1 = gc(7), th_2 = gc(8), th_3 = gc(9);

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

    if (i==0)
        end_pos = r_0;
    else if (i==1)
        end_pos = r_0 + Rot_0*r_1;
    else if (i==2)
        end_pos = r_0 + Rot_0*(r_1 + Rot_1*r_2);
    else if (i==3)
        end_pos = r_0 + Rot_0*(r_1 + Rot_1*(r_2 + Rot_2*r_3));
    else if (i==4)
        end_pos = r_0 + Rot_0*(r_1 + Rot_1*(r_2 + Rot_2*(r_3 + Rot_3 * r_4)));


    return end_pos; /// replace this
}


/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Eigen::Vector3d pos_root, pos_joint_1, pos_joint_2, pos_joint_3, pos_joint_4, linear_v;
    Eigen::Vector3d r_0, r_1, r_2, r_3, v_0, v_1, v_2, v_3, v_4;

    pos_root = getEndEffectorPosition(gc, 0);
    pos_joint_1 = getEndEffectorPosition(gc, 1);
    pos_joint_2 = getEndEffectorPosition(gc, 2);
    pos_joint_3 = getEndEffectorPosition(gc, 3);
    pos_joint_4 = getEndEffectorPosition(gc, 4);

    r_0 = pos_joint_4 - pos_root;
    r_1 = pos_joint_4 - pos_joint_1;
    r_2 = pos_joint_4 - pos_joint_2;
    r_3 = pos_joint_4 - pos_joint_3;



    return linear_v; /// replace this
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Eigen::Vector3d pos_joint_1, pos_joint_2, pos_joint_3, pos_joint_4, angular_v;

    pos_joint_1 = getEndEffectorPosition(gc, 1);
    pos_joint_2 = getEndEffectorPosition(gc, 2);
    pos_joint_3 = getEndEffectorPosition(gc, 3);
    pos_joint_4 = getEndEffectorPosition(gc, 4);

  return angular_v; /// replace this
}