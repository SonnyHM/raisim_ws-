#pragma once
#include <iostream>
#include <memory>
#include <Eigen/Core>
#include "Eigen/Geometry"
#include <cmath>
#include <list>
//

Eigen::Matrix3d skew(const Eigen::Vector3d &vt){
    Eigen::Matrix3d cross_matrix;
    cross_matrix <<0, -vt(2), vt(1),
            vt(2), 0, -vt(0),
            -vt(1),vt(0),0;
    return cross_matrix;
}

namespace aliengo{
    class Body{
    public:
        struct Joint{
            /// variables
            Eigen::Vector3d jointPosition_W;
            Eigen::Vector3d jointAxis_W;
            Eigen::Matrix3d rotation_W;


            /// robot definition
            Eigen::Vector3d jointPosition_P;
            Eigen::Vector3d jointAxis_P;
            Eigen::Matrix3d jointRotations_P;
            enum class Type{
                FIXED,
                FLOATING,
                REVOLUTE,
                PRISMATIC
            } type = Type::REVOLUTE;
        };

        Eigen::Vector3d COM_pose;
        Eigen::Matrix3d Inertia_W;
        Eigen::VectorXd motion_S;
        double Mass;

        Body(const Joint::Type type,
             const Eigen::Vector3d& jointAxis_P,
             const Eigen::Vector3d& jointPosition_P,
             const Eigen::Matrix3d& jointRotations_P,
             const Eigen::Vector3d& bodyframe_P,
             const Eigen::Matrix3d& Inertia_P,
             const double link_mass){

          joint_.jointAxis_P = jointAxis_P;
          joint_.jointPosition_P = jointPosition_P;
          joint_.jointRotations_P = jointRotations_P;
          joint_.type = type;
          COM_pose = bodyframe_P;
          Inertia_W= Inertia_P;
          Mass = link_mass;
          motion_S.resize(6);
        }

        void computeKinematicsDownTheTree(const Eigen::VectorXd& gc){

          std::vector<Eigen::VectorXd> gc_list;
          Eigen::VectorXd gc_j;
          if (joint_.type== Joint::Type::FLOATING){
            ///joint position
              Eigen::Vector3d p_pose(gc(0), gc(1), gc(2));
              Eigen::Quaterniond parent_q(gc(3),gc(4),gc(5),gc(6));

              joint_.rotation_W = parent_q.normalized().toRotationMatrix();
              joint_.jointPosition_W = p_pose;
            for (int i=0; i <children_.size(); i++){
              gc_j=gc.segment(7+3*i,3);
              gc_list.push_back(gc_j);
            }
            ///velocity
          }

          else if(joint_.type== Joint::Type::REVOLUTE){
            ///find joint position
            Eigen::Matrix3d Rotation;
            if (joint_.jointAxis_P == Eigen::Vector3d (1,0,0)){
              Rotation << 1, 0, 0,
                      0, cos(gc(0)), -sin(gc(0)),
                      0, sin(gc(0)), cos(gc(0));
            }
            else if (joint_.jointAxis_P == Eigen::Vector3d (0,1,0)){
              Rotation<< cos(gc(0)), 0, sin(gc(0)),
                      0, 1, 0,
                      -sin(gc(0)), 0, cos(gc(0));
            }
            else {
              Rotation<< cos(gc(0)), -sin(gc(0)), 0,
                      sin(gc(0)), cos(gc(0)), 0,
                      0,0,1;
            }
            joint_.jointPosition_W = parent_->joint_.jointPosition_W + parent_->joint_.rotation_W*joint_.jointPosition_P;
            joint_.rotation_W = parent_->joint_.rotation_W * joint_.jointRotations_P* Rotation;
            joint_.jointAxis_W = (joint_.rotation_W*joint_.jointAxis_P).normalized();
            motion_S << 0,0,0, joint_.jointAxis_W;
              gc_j=gc.segment(1,gc.size()-1);
              gc_list.push_back(gc_j);
          }

          else if(joint_.type== Joint::Type::PRISMATIC){
              Eigen::Vector3d prismatic_pos;
              if (joint_.jointAxis_P == Eigen::Vector3d (1,0,0)){
                  prismatic_pos << gc(0), 0, 0;
                  joint_.jointPosition_W = parent_->joint_.jointPosition_W + parent_->joint_.rotation_W*(prismatic_pos + joint_.jointPosition_P);
                  joint_.rotation_W = parent_->joint_.rotation_W * joint_.jointRotations_P;
              }
              else if (joint_.jointAxis_P == Eigen::Vector3d (0,1,0))
              {
                  prismatic_pos << 0, gc(0), 0;
                  joint_.jointPosition_W = parent_->joint_.jointPosition_W + parent_->joint_.rotation_W*( prismatic_pos+ joint_.jointPosition_P);
                  joint_.rotation_W = parent_->joint_.rotation_W * joint_.jointRotations_P;
              }
              joint_.jointAxis_W=(parent_->joint_.rotation_W*joint_.jointAxis_P).normalized();
              motion_S << joint_.jointAxis_W,0,0,0;
              gc_j=gc.segment(1,gc.size()-1);
              gc_list.push_back(gc_j);

          }

          else if(joint_.type== Joint::Type::FIXED){
            joint_.jointPosition_W = parent_->joint_.jointPosition_W+ parent_->joint_.rotation_W*joint_.jointPosition_P;
            joint_.rotation_W = parent_->joint_.rotation_W* joint_.jointRotations_P;
          }

          COM_pose= joint_.jointPosition_W + joint_.rotation_W*COM_pose;
          Inertia_W= joint_.rotation_W*Inertia_W*joint_.rotation_W.transpose();

          int i = 0;
          for (auto &child: children_) {
              child->parent_ = this;
              child->computeKinematicsDownTheTree(gc_list[i]);
              i = i + 1;
          }
        }

        void setChildren(const std::vector<Body*> children){
          children_ = children;
        }

        Joint joint_;
        std::vector<Body*> children_;
        Body* parent_ ;
    };

    class composite_body{
    public:
        std::vector<Body*> joint_order;
        Eigen::Vector3d COM;
        Eigen::Matrix3d Composite_I;
        double Composite_mass;
        std::vector<Eigen::Matrix3d> list_I;
        std::vector<Eigen::Vector3d> list_COM;
        std::vector<double> list_Mass;


        composite_body(const std::vector<Body*> to_joint_order){
            joint_order = to_joint_order;
        }

        void compute_composite_body() {
            Eigen::Vector3d r_1, r_2;
            Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
            double m_1,m_2;

            for (int i = joint_order.size()-1; i > 0;i--) {
                if (i==joint_order.size()-1){
                    m_1= joint_order[i-1]->Mass;
                    m_2 = joint_order[i]->Mass;
                    COM = (m_1*joint_order[i-1]->COM_pose+m_2*joint_order[i]->COM_pose)/(m_1+m_2);
                    r_1= joint_order[i-1]->COM_pose - COM;
                    r_2= joint_order[i]->COM_pose - COM;
                    Composite_I=joint_order[i-1]->Inertia_W + joint_order[i]->Inertia_W - m_1*I*skew(r_1)*skew(r_1)- m_2*I*skew(r_2)*skew(r_2);
                }
                else{
                    m_1= joint_order[i-1]->Mass;
                    m_2 = list_Mass[0];
                    COM = (m_1*joint_order[i-1]->COM_pose+m_2*list_COM[0])/(m_1+m_2);
                    r_1= joint_order[i-1]->COM_pose - COM;
                    r_2= list_COM[0] - COM;
                    Composite_I=joint_order[i-1]->Inertia_W + list_I[0] - m_1*skew(r_1)*skew(r_1)- m_2*skew(r_2)*skew(r_2);
                }

                Composite_mass=m_1+m_2;
                list_I.insert(list_I.begin(),Composite_I);
                list_COM.insert(list_COM.begin(),COM);
                list_Mass.insert(list_Mass.begin(),Composite_mass);

            }
        }
    };
}

/// do not change the name of the method

inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {
    Eigen::MatrixXd mass_matrix(18,18), S_Mass_i(6,6);
    Eigen::Vector3d pos_root(gc(0),gc(1),gc(2)),com_root,com_root_prev;
    Eigen::Vector3d r_0, r_1, r_2, r_3, v_0, w_0, w_1, w_2, w_3;
    Eigen::Matrix3d Rot_I, Inertia_IMU,Inertia_0, Inertia_1, Inertia_2, Inertia_3, Inertia_4,root_Inertia, Inertia_5, Inertia_6, Inertia_7,Inertia_8,Inertia_9,Inertia_10,Inertia_11,Inertia_12,Inertia_13,Inertia_14,Inertia_15,Inertia_16;
    Eigen::Quaterniond parent_q(gc(3),gc(4),gc(5),gc(6));
    std::vector<aliengo::Body> leg1,leg2,leg3,leg4;
    double base_m, m_IMU, m_0,m_1, m_2, m_3, m_4;
    std::vector<Eigen::VectorXd> trunk_m_Matrix;
    std::vector<double> S_Mass;
    std::vector<aliengo::composite_body*> articulated;
    mass_matrix.setZero();
    Eigen::MatrixXd Last_S_Mass(6,6);
    Eigen::MatrixXd trunk_motion_S = Eigen::MatrixXd::Identity(6, 6);

    m_IMU = 0.001;
    m_0 = 9.041;
    m_1= 1.993;
    m_2=0.639;
    m_3=0.207;
    m_4= 0.06;

    Rot_I << 1,0,0,
            0,1,0,
            0,0,1;

    Inertia_IMU<< 0.0001,0,0,
                  0,0.000001,0,
                  0,0,0.0001;

    Inertia_0 <<0.033260231,-0.000451628,  0.000487603,
            -0.000451628, 0.16117211, 4.8356e-05,
            0.000487603 ,4.8356e-05,0.17460442;

    Inertia_1 <<0.002903894, 7.185e-05,-1.262e-06,
            7.185e-05,0.004907517, 1.75e-06,
            -1.262e-06, 1.75e-06,0.005586944;

    Inertia_2 <<0.005666803,-3.597e-06,  0.000491446,
            -3.597e-06, 0.005847229, -1.0086e-05,
            0.000491446 ,-1.0086e-05, 0.000369811;

    Inertia_3 <<0.006341369,-3e-09,  -8.7951e-05,
            -3e-09, 0.006355157, -1.336e-06,
            -8.7951e-05 ,-1.336e-06,3.9188e-05;

    Inertia_4 <<1.6854e-05,0, 0,
            0, 1.6854e-05, 0,
            0,0,1.6854e-05;

    Inertia_5 <<0.002903894, -7.185e-05,-1.262e-06,
            -7.185e-05,0.004907517, -1.75e-06,
            -1.262e-06, -1.75e-06,0.005586944;

    Inertia_6 <<0.005666803,3.597e-06,  0.000491446,
            3.597e-06, 0.005847229, 1.0086e-05,
            0.000491446 ,1.0086e-05, 0.000369811;

    Inertia_7 <<0.006341369,-3e-09,  -8.7951e-05,
            -3e-09, 0.006355157, -1.336e-06,
            -8.7951e-05 ,-1.336e-06,3.9188e-05;

    Inertia_8 <<1.6854e-05,0, 0,
            0, 1.6854e-05, 0,
            0,0,1.6854e-05;

    Inertia_9 <<0.002903894, -7.185e-05,1.262e-06,
            -7.185e-05,0.004907517, 1.75e-06,
            1.262e-06, 1.75e-06,0.005586944;

    Inertia_10 <<0.005666803,-3.597e-06,  0.000491446,
            -3.597e-06, 0.005847229, -1.0086e-05,
            0.000491446 ,-1.0086e-05, 0.000369811;

    Inertia_11 <<0.006341369,-3e-09,  -8.7951e-05,
            -3e-09, 0.006355157, -1.336e-06,
            -8.7951e-05 ,-1.336e-06,3.9188e-05;

    Inertia_12 <<1.6854e-05,0, 0,
            0, 1.6854e-05, 0,
            0,0,1.6854e-05;

    Inertia_13 <<0.002903894, 7.185e-05,1.262e-06,
            7.185e-05,0.004907517, -1.75e-06,
            1.262e-06, -1.75e-06,0.005586944;

    Inertia_14 <<0.005666803,3.597e-06,  0.000491446,
            3.597e-06, 0.005847229, 1.0086e-05,
            0.000491446 ,1.0086e-05, 0.000369811;

    Inertia_15 <<0.006341369,-3e-09,  -8.7951e-05,
            -3e-09, 0.006355157, -1.336e-06,
            -8.7951e-05 ,-1.336e-06,3.9188e-05;

    Inertia_16 <<1.6854e-05,0, 0,
            0, 1.6854e-05, 0,
            0,0,1.6854e-05;

    aliengo::Body* IMU= new aliengo::Body (aliengo::Body::Joint::Type::FIXED,
                                             {0,0,0},
                                             {0,0,0},
                                             Rot_I,
                                             {0, 0, 0 },
                                             Inertia_IMU,
                                             m_IMU);

    aliengo::Body* trunk_1= new aliengo::Body (aliengo::Body::Joint::Type::FLOATING,
                                             {0,0,0},
                                             {0,0,0},
                                             Rot_I,
                                             {0.008465, 0.004045, -0.000763},
                                             Inertia_0,
                                             m_0);

    aliengo::Body* trunk= new aliengo::Body (aliengo::Body::Joint::Type::FLOATING,
                                             {0,0,0},
                                             {0,0,0},
                                             Rot_I,
                                             {0.008465, 0.004045, -0.000763},
                                             Inertia_0,
                                             m_0);


    aliengo::Body* FR_hip = new aliengo::Body(aliengo::Body::Joint::Type::REVOLUTE,
                                              {1,0,0},
                                              {0.2399,-0.051, 0},
                                              Rot_I,
                                              {-0.022191, -0.015144, -0.000015},
                                              Inertia_1,
                                              m_1);

    aliengo::Body* FR_thigh = new aliengo::Body(aliengo::Body::Joint::Type::REVOLUTE,
                                                {0,1,0},
                                                {0, -0.083, 0},
                                                Rot_I,
                                                {-0.005607, 0.003877, -0.048199}
                                                ,Inertia_2,
                                                m_2);

    aliengo::Body* FR_calf = new aliengo::Body(aliengo::Body::Joint::Type::REVOLUTE,
                                               {0,1,0},
                                               {0, 0, -0.25},
                                               Rot_I,
                                               {0.002781, 6.3e-05, -0.142518},
                                               Inertia_3,
                                               m_3);

    aliengo::Body* FR_foot_fixed = new aliengo::Body(aliengo::Body::Joint::Type::FIXED,
                                                     {0,0,0},
                                                     {0, 0, -0.25},
                                                     Rot_I,
                                                     {0,0,0},
                                                     Inertia_4,
                                                     m_4);

    aliengo::Body* FL_hip = new aliengo::Body(aliengo::Body::Joint::Type::REVOLUTE,
                                              {1,0,0},
                                              {0.2399,0.051, 0},
                                              Rot_I,
                                              {-0.022191, 0.015144, -1.5e-05},
                                              Inertia_5,
                                              m_1);

    aliengo::Body* FL_thigh = new aliengo::Body(aliengo::Body::Joint::Type::PRISMATIC,
                                                {0,1,0},
                                                {0, 0.083, 0},
                                                Rot_I,
                                                {-0.005607, -0.003877, -0.048199}
                                                ,Inertia_6,
                                                m_2);

    aliengo::Body* FL_calf = new aliengo::Body(aliengo::Body::Joint::Type::REVOLUTE,
                                               {0,1,0},
                                               {0, 0, -0.25},
                                               Rot_I,
                                               {0.002781, 6.3e-05, -0.142518},
                                               Inertia_7,
                                               m_3);

    aliengo::Body* FL_foot_fixed = new aliengo::Body(aliengo::Body::Joint::Type::FIXED,
                                                     {0,0,0},
                                                     {0, 0, -0.25},
                                                     Rot_I,
                                                     {0,0,0},
                                                     Inertia_8,
                                                     m_4);

    aliengo::Body* RR_hip = new aliengo::Body(aliengo::Body::Joint::Type::REVOLUTE,
                                              {1,0,0},
                                              {-0.2399,-0.051, 0},
                                              Rot_I,
                                              {0.022191, -0.015144, -0.000015},
                                              Inertia_9,
                                              m_1);

    aliengo::Body* RR_thigh = new aliengo::Body(aliengo::Body::Joint::Type::REVOLUTE,
                                                {0,1,0},
                                                {0, -0.083, 0},
                                                Rot_I,
                                                {-0.005607, 0.003877, -0.048199}
                                                ,Inertia_10,
                                                m_2);

    aliengo::Body* RR_calf = new aliengo::Body(aliengo::Body::Joint::Type::REVOLUTE,
                                               {0,1,0},
                                               {0, 0, -0.25},
                                               Rot_I,
                                               {0.002781, 6.3e-05, -0.142518},
                                               Inertia_11,
                                               m_3);

    aliengo::Body* RR_foot_fixed = new aliengo::Body(aliengo::Body::Joint::Type::FIXED,
                                                     {0,0,0},
                                                     {0, 0, -0.25},
                                                     Rot_I,
                                                     {0,0,0},
                                                     Inertia_12,
                                                     m_4);

    aliengo::Body* RL_hip = new aliengo::Body(aliengo::Body::Joint::Type::PRISMATIC,
                                              {1,0,0},
                                              {-0.2399,0.051, 0},
                                              Rot_I,
                                              {0.022191, 0.015144, -1.5e-05},
                                              Inertia_13,
                                              m_1);

    aliengo::Body* RL_thigh = new aliengo::Body(aliengo::Body::Joint::Type::PRISMATIC,
                                                {0,1,0},
                                                {0, 0.083, 0},
                                                Rot_I,
                                                {-0.005607, -0.003877, -0.048199}
                                                ,Inertia_14,
                                                m_2);

    aliengo::Body* RL_calf = new aliengo::Body(aliengo::Body::Joint::Type::PRISMATIC,
                                               {0,1,0},
                                               {0, 0, -0.25},
                                               Rot_I,
                                               {0.002781, 6.3e-05, -0.142518},
                                               Inertia_15,
                                               m_3);

    aliengo::Body* RL_foot_fixed = new aliengo::Body(aliengo::Body::Joint::Type::FIXED,
                                                     {0,0,0},
                                                     {0, 0, -0.25},
                                                     Rot_I,
                                                     {0,0,0},
                                                     Inertia_16,
                                                     m_4);

    FR_calf->setChildren({FR_foot_fixed});
    FR_thigh->setChildren({FR_calf});
    FR_hip->setChildren({FR_thigh});
    FL_calf->setChildren({FL_foot_fixed});
    FL_thigh->setChildren({FL_calf});
    FL_hip->setChildren({FL_thigh});
    RR_calf->setChildren({RR_foot_fixed});
    RR_thigh->setChildren({RR_calf});
    RR_hip->setChildren({RR_thigh});
    RL_calf->setChildren({RL_foot_fixed});
    RL_thigh->setChildren({RL_calf});
    RL_hip->setChildren({RL_thigh});
    trunk->setChildren({FR_hip, FL_hip, RR_hip, RL_hip});
    trunk->computeKinematicsDownTheTree(gc);


    trunk_1->setChildren({IMU});
    trunk_1->computeKinematicsDownTheTree(gc);


    aliengo::composite_body* base = new aliengo::composite_body({IMU, trunk_1});

    base->compute_composite_body();

    trunk -> COM_pose= base->COM;
    trunk -> Mass = base ->Composite_mass;
    trunk -> Inertia_W = base-> Composite_I;

    aliengo::composite_body* FR_leg = new aliengo::composite_body({trunk,FR_hip,FR_thigh,FR_calf,FR_foot_fixed});
    aliengo::composite_body* FL_leg = new aliengo::composite_body({trunk,FL_hip,FL_thigh,FL_calf,FL_foot_fixed});
    aliengo::composite_body* RR_leg = new aliengo::composite_body({trunk,RR_hip,RR_thigh,RR_calf,RR_foot_fixed});
    aliengo::composite_body* RL_leg = new aliengo::composite_body({trunk,RL_hip,RL_thigh,RL_calf,RL_foot_fixed});

    FR_leg -> compute_composite_body();
    FL_leg -> compute_composite_body();
    RL_leg -> compute_composite_body();
    RR_leg -> compute_composite_body();


    S_Mass_i.setZero();
    mass_matrix.setZero();
    Last_S_Mass.setZero();
    articulated = {FR_leg,FL_leg,RR_leg,RL_leg};

    com_root= trunk->COM_pose;
    root_Inertia = trunk->Inertia_W;
    m_0 = trunk->Mass;

    for(int k=3; k>-1;k--) {
        for (int j = 3; j > 0; j--) {
            for (int i = j; i > -1; i--) {
                Eigen::Vector3d r_ac = articulated[k]->list_COM[j] - articulated[k]->joint_order[j]->joint_.jointPosition_W;
                Eigen::MatrixXd R = Eigen::MatrixXd::Identity(6, 6);
                R.block(0, 3, 3, 3) = -skew(articulated[k]->joint_order[j]->joint_.jointPosition_W-articulated[k]->joint_order[i]->joint_.jointPosition_W);
                S_Mass_i.block(0, 0, 3, 3) = articulated[k]->list_Mass[j] * Eigen::Matrix3d::Identity();
                S_Mass_i.block(3, 0, 3, 3) = articulated[k]->list_Mass[j] * skew(r_ac);
                S_Mass_i.block(0, 3, 3, 3) = -articulated[k]->list_Mass[j] * skew(r_ac);
                S_Mass_i.block(3, 3, 3, 3) = articulated[k]->list_I[j] - articulated[k]->list_Mass[j] * skew(r_ac) * skew(r_ac);
                S_Mass_i = S_Mass_i * R;
                if(i==0){
                    mass_matrix.block(0,5+3*k+j,6,1)= (articulated[k]->joint_order[j]->motion_S.transpose()*S_Mass_i).transpose();
                    mass_matrix.block(5+3*k+j,0,1,6)= articulated[k]->joint_order[j]->motion_S.transpose()*S_Mass_i;

                }
                else{
                    mass_matrix.block(5+3*k+j,5+3*k+i,1,1) = articulated[k]->joint_order[j]->motion_S.transpose()*S_Mass_i*articulated[k]->joint_order[i]->motion_S;
                    mass_matrix.block(5+3*k+i,5+3*k+j,1,1) = articulated[k]->joint_order[j]->motion_S.transpose()*S_Mass_i*articulated[k]->joint_order[i]->motion_S;
                }
            }
        }
        com_root_prev = com_root;
        com_root= (m_0*com_root_prev+articulated[k]->list_Mass[1]*articulated[k]->list_COM[1])/(m_0 + articulated[k]->list_Mass[1]);
        root_Inertia=root_Inertia+articulated[k]->list_I[1]-m_0*skew(com_root_prev-com_root)*skew(com_root_prev-com_root)-articulated[k]->list_Mass[1]*skew(articulated[k]->list_COM[1]-com_root)*skew(articulated[k]->list_COM[1]-com_root);
        m_0= m_0 + articulated[k]->list_Mass[1];
    }

    Last_S_Mass.block(0,0,3,3) = m_0*Eigen::Matrix3d::Identity();
    Last_S_Mass.block(3,0,3,3) = m_0*skew(com_root-trunk->joint_.jointPosition_W);
    Last_S_Mass.block(0,3,3,3) = -m_0*skew(com_root-trunk->joint_.jointPosition_W);
    Last_S_Mass.block(3,3,3,3) = root_Inertia-m_0*skew(com_root-trunk->joint_.jointPosition_W)*skew(com_root-trunk->joint_.jointPosition_W);

    mass_matrix.block(0,0,6,6) = Last_S_Mass;




  return mass_matrix;
}