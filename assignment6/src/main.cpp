#include <iostream>
#include <math.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any libigl headers here ***/
#include <igl/readTGF.h>
#include <igl/readDMAT.h>
#include <igl/directed_edge_parents.h>
#include <igl/column_to_quats.h>
#include <igl/deform_skeleton.h>
#include <igl/forward_kinematics.h>
/*** insert any libhedra headers here ***/
#include <hedra/line_cylinders.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

typedef
  std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond>>
  RotationList;

// for animation
RotationList pose;
double anim_t = 1.0;
double anim_t_dir = -0.03;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// skeleton vertices by 3 list of vertex positions
Eigen::MatrixXd C;
// skeleton edges by 2 list of edge indices
Eigen::MatrixXi E;
// list of parents
Eigen::VectorXi P;
// transformation matrices
Eigen::MatrixXd transM;
// transformation quaternions
Eigen::MatrixXd transQ;

bool mesh_loaded = false;
// load mesh
bool load_mesh(Viewer& viewer, string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
  igl::readOFF(filename, V, F);
  viewer.data().clear();
  viewer.data().set_mesh(V, F);
  viewer.data().compute_normals();
  viewer.core().align_camera_center(V, F);

  return true;
}

// show mesh in the viewer
void show_mesh(Viewer& viewer, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.data().compute_normals();
    viewer.core().align_camera_center(V, F);

    return ;
}

// load and show skeleton in the viewer
bool load_skeleton(Viewer& viewer, string filename, Eigen::MatrixXd& C, Eigen::MatrixXi& E) {
    // read .tgf
    igl::readTGF(filename, C, E);

    // compute line cylinders given C and E
    // compute P1 and p2
    Eigen::MatrixXd P1, P2;
    P1.resize(E.rows(), 3);
    P2.resize(E.rows(), 3);
    for (int i = 0; i < E.rows(); i ++) {
        P1.row(i) << C.row(E(i, 0));
        P2.row(i) << C.row(E(i, 1));
    }

    // compute ptColors
    Eigen::MatrixXd ptColors;
    ptColors.resize(E.rows(), 3);
    for (int i = 0; i < E.rows(); i ++) {
        ptColors.row(i) << 1.0, 0.0, 0;
    }
    // compute cyndColors
    Eigen::MatrixXd cyndColors;
    cyndColors.resize(E.rows(), 3);
    for (int i = 0; i < E.rows(); i ++) {
        cyndColors.row(i) << 0.0, 1.0, 0;
    }

    // compute line cylinders
    Eigen::MatrixXd cyndV, cyndC;
    Eigen::MatrixXi cyndF;
    hedra::line_cylinders(P1, P2, 0.01, cyndColors, 10, cyndV, cyndF, cyndC);

    // add data to viewer
    viewer.data().clear();
    viewer.data().set_points(C, ptColors);
    viewer.data().set_mesh(cyndV, cyndF);
    viewer.data().set_colors(cyndC);
    viewer.data().compute_normals();
    viewer.core().align_camera_center(cyndV, cyndF);

    return true;
}

// void show_skeleton()

bool load_matrot(string filename, Eigen::MatrixXd& transM) {
    igl::readDMAT(filename, transM);

    return true;
}

bool load_handles() {
    return true;
}

bool load_quat(string filename, Eigen::MatrixXd& transQ) {
    igl::readDMAT(filename, transQ);

    return true;
}

void compute_absolute_transmat() {
    // TODO

    return ;
}

void compute_absolute_transquat() {
    

    return ;
}

bool pre_draw (igl::opengl::glfw::Viewer & viewer) {
    using namespace Eigen;

    if(viewer.core().is_animating) {
        // Interpolate pose and identity
        RotationList anim_pose(pose.size());
        for(int e = 0; e < pose.size(); e++) {
            anim_pose[e] = pose[e].slerp(anim_t,Quaterniond::Identity());
        }

        // Propagate relative rotations via FK to retrieve absolute transformations
        RotationList vQ;
        vector<Vector3d> vT;
        igl::forward_kinematics(C, E, P, anim_pose, vQ, vT);
        const int dim = C.cols();
        MatrixXd T(E.rows()*(dim+1),dim);
        for(int e = 0; e < E.rows(); e++) {
            Affine3d a = Affine3d::Identity();
            a.translate(vT[e]);
            a.rotate(vQ[e]);
            T.block(e*(dim+1),0,dim+1,dim) =
            a.matrix().transpose().block(0,0,dim+1,dim);
        }

        // Also deform skeleton edges
        MatrixXd CT;
        MatrixXi BET;
        igl::deform_skeleton(C, E, T, CT, BET);

        // compute cylinders
        // compute P1 and p2
        Eigen::MatrixXd P1, P2;
        P1.resize(E.rows(), 3);
        P2.resize(E.rows(), 3);
        for (int i = 0; i < E.rows(); i ++) {
            P1.row(i) << CT.row(BET(i, 0));
            P2.row(i) << CT.row(BET(i, 1));
        }

        // compute ptColors
        Eigen::MatrixXd ptColors;
        ptColors.resize(E.rows(), 3);
        for (int i = 0; i < E.rows(); i ++) {
            ptColors.row(i) << 1.0, 0.0, 0;
        }
        // compute cyndColors
        Eigen::MatrixXd cyndColors;
        cyndColors.resize(E.rows(), 3);
        for (int i = 0; i < E.rows(); i ++) {
            cyndColors.row(i) << 0.0, 1.0, 0;
        }

        // compute line cylinders
        Eigen::MatrixXd cyndV, cyndC;
        Eigen::MatrixXi cyndF;
        hedra::line_cylinders(P1, P2, 0.01, cyndColors, 10, cyndV, cyndF, cyndC);

        viewer.data().set_points(CT, ptColors);
        viewer.data().set_mesh(cyndV, cyndF);
        viewer.data().set_colors(cyndC);
        viewer.data().compute_normals();
        anim_t += anim_t_dir;
        anim_t_dir *= (anim_t>=1.0 || anim_t<=0.0?-1.0:1.0);
    }

    return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    case ' ':
      viewer.core().is_animating = !viewer.core().is_animating;
      break;
  }

  return true;
}

int main(int argc, char *argv[]) {
    // Show the mesh
    Viewer viewer;
    
    std::string filename;
    if (argc == 2) {
        filename = std::string(argv[1]); // Mesh provided as command line argument
    }
    else {
        filename = std::string("../data/hand/hand.tgf"); // Default mesh
    }
	
    // mesh_loaded = load_mesh(viewer, filename, V, F);
    // show_mesh(viewer, V, F);
    load_skeleton(viewer, filename, C, E);
    igl::directed_edge_parents(E, P);
    load_quat("../data/hand/hand-pose_quat.dmat", transQ);
    igl::column_to_quats(transQ, pose);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &key_down;
    viewer.core().is_animating = false;
    viewer.core().animation_max_fps = 30.;
    
    viewer.launch();
}