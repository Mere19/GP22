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
#include <igl/readOBJ.h>
#include <igl/directed_edge_parents.h>
#include <igl/column_to_quats.h>
#include <igl/deform_skeleton.h>
#include <igl/forward_kinematics.h>
/*** insert any libhedra headers here ***/
#include <hedra/line_cylinders.h>
/*** insert any self-defined headers here ***/
#include "gui.h"
#include "control.h"
#include "handles.h"
#include "compute.h"
// #include "utils.h"

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// for handles
int handle_id = 0;
Eigen::VectorXi ref_handle_vertices;
Eigen::VectorXi ref_free_vertices;
Eigen::VectorXi ref_root_vertices;
Eigen::VectorXi ref_non_root_vertices;
Eigen::MatrixXd ref_root_handle_positions;
Eigen::VectorXi handle_vertices;
Eigen::VectorXi free_vertices;
Eigen::VectorXi root_vertices;
Eigen::VectorXi non_root_vertices;
Eigen::MatrixXd root_handle_positions;

// for skinning weight function
Eigen::MatrixXd W, FW, fQG;
Eigen::SparseMatrix<double> Lw, G, D, lAff, lAfc, pAff, pAfc;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> laplace_solver;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> poisson_solver;

// for animation
RotationList pose;
int frame = 0;
double anim_t = 1.0;
double anim_t_dir = -0.03;
Eigen::MatrixXd Rot, Trans;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// skeleton vertices by 3 list of vertex positions
Eigen::MatrixXd C;
// skeleton edges by 2 list of edge indices
Eigen::MatrixXi E;
// #V x1 reference handle ID, -1 = not handle vertices
Eigen::VectorXi refH;
// #V x1 handle ID, -1 = not handle vertices
Eigen::VectorXi H;
// list of parents
Eigen::VectorXi P;
// transformation matrices
Eigen::MatrixXd transM;
// transformation quaternions
Eigen::MatrixXd transQ;
// handle IDs
Eigen::MatrixXd handleID;

// flags
bool mesh_loaded = false;
bool skeleton_loaded = false;
bool matrot = false;
bool visualize_skeleton = false;
bool visualize_mesh = false;

// load mesh
bool load_mesh (string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    igl::readOFF(filename, V, F);

    return true;
}

bool load_skeleton (string filename, Eigen::MatrixXd& C, Eigen::MatrixXi& E) {
    // read .tgf
    igl::readTGF(filename, C, E);

    return true;
}

bool load_matrot (string filename, Eigen::MatrixXd& transM) {
    igl::readDMAT(filename, transM);

    return true;
}

bool load_handles (string filename, Eigen::VectorXi& H) {
    igl::readDMAT(filename, H);
    
    return true;
}

bool load_quat (string filename, Eigen::MatrixXd& transQ) {
    igl::readDMAT(filename, transQ);

    return true;
}

void compute_absolute_transmat(Eigen::MatrixXd& transM, int& frame, Eigen::MatrixXd& T) {
    using namespace Eigen;

    const int dim = C.cols();
    T.resize(E.rows() * (dim + 1), dim);

    // propagate relative rotations via FK to retrieve absolute transformations
    vector<Matrix3d> vQ;
    vector<Vector3d> vT;
    MatrixXd dQ = transM.block(frame * E.rows() * dim, 0, E.rows() * dim, dim);
    forward_kinematics(C, E, P, dQ, vQ, vT);
    for (int e = 0; e < E.rows(); e ++) {
        Affine3d a = Affine3d::Identity();
        a.translate(vT[e]);
        a.rotate(vQ[e]);        // support dim * dim rotation matrix as well
        T.block(e * (dim + 1), 0, dim + 1, dim) = a.matrix().transpose().block(0, 0, dim + 1, dim);
    }

    return;
}

// void compute_absolute_transmat(Eigen::MatrixXd& transM, int& frame, Eigen::MatrixXd& T) {
//     using namespace Eigen;

//     const int dim = C.cols();
//     T.resize(E.rows() * (dim + 1), dim);

//     // propagate relative rotations via FK to retrieve absolute transformations
//     MatrixXd dQ = transM.block(frame * E.rows() * dim, 0, E.rows() * dim, dim);
//     vector<Affine3d> AList;
//     forward_kinematics__(C, E, P, dQ, AList);
//     for (int e = 0; e < E.rows(); e ++) {
//         T.block(e * (dim + 1), 0, dim + 1, dim) = AList[e].matrix().transpose().block(0, 0, dim + 1, dim);
//     }

//     return;
// }

void compute_absolute_transquat (RotationList& pose, double& anim_t, Eigen::MatrixXd& T) {
    using namespace Eigen;

    const int dim = C.cols();
    T.resize(E.rows() * (dim + 1), dim);

    // interpolate pose and identity
    RotationList anim_pose(pose.size());
    for(int e = 0; e < pose.size(); e++) {
        anim_pose[e] = pose[e].slerp(anim_t, Quaterniond::Identity());
    }

    // propagate relative rotations via FK to retrieve absolute transformations
    RotationList vQ;
    vector<Vector3d> vT;
    igl::forward_kinematics(C, E, P, anim_pose, vQ, vT);
    for (int e = 0; e < E.rows(); e ++) {
        Affine3d a = Affine3d::Identity();
        a.translate(vT[e]);
        a.rotate(vQ[e]);        // support dim * dim rotation matrix as well
        T.block(e * (dim + 1), 0, dim + 1, dim) = a.matrix().transpose().block(0, 0, dim + 1, dim);
    }

    return ;
}

void compute_absolute_RTquat (RotationList& pose, double& anim_t, RotationList& vQ, vector<Eigen::Vector3d>& vT) {
    using namespace Eigen;

    const int dim = C.cols();

    // interpolate pose and identity
    RotationList anim_pose(pose.size());
    for(int e = 0; e < pose.size(); e++) {
        anim_pose[e] = pose[e].slerp(anim_t, Quaterniond::Identity());
    }

    // propagate relative rotations via FK to retrieve absolute transformations
    igl::forward_kinematics(C, E, P, anim_pose, vQ, vT);
    // for (int e = 0; e < E.rows(); e ++) {
    //     Affine3d a = Affine3d::Identity();
    //     a.translate(vT[e]);
    //     a.rotate(vQ[e]);        // support dim * dim rotation matrix as well
    //     T.block(e * (dim + 1), 0, dim + 1, dim) = a.matrix().transpose().block(0, 0, dim + 1, dim);
    // }

    return ;
}

bool pre_draw (igl::opengl::glfw::Viewer & viewer) {
    using namespace Eigen;

    RotationList vQ;
    vector<Vector3d> vT;

    if(viewer.core().is_animating) {
        // compute absolute transformations
        MatrixXd T;
        // if (matrot) {
            compute_absolute_transmat(transM, frame, T);
        // } else {
            // compute_absolute_transquat(pose, anim_t, T);
            compute_absolute_RTquat(pose, anim_t, vQ, vT);
        // }

        // deform skeleton
        MatrixXd CT;
        MatrixXi BET;
        // igl::deform_skeleton(C, E, T, CT, BET);

        // per-vertex linear blending skinning
        MatrixXd VT;
        // compute_per_vertex_linear_blending_skinning(V, W, T, VT);

        // dual quaternion skinning
        // compute_dual_quaternion_skinning(V, W, vQ, vT, VT);

        // per-face linear blending skinning
        RotationList fQ;
        compute_per_face_transformation(FW, vQ, fQ);
        compute_face_deformation_gradient(fQ, fQG);
        compute_poisson_stitching(V, G, D, fQG, poisson_solver, pAfc, ref_root_vertices, ref_non_root_vertices, ref_root_handle_positions, VT);

        // display deformed skeleton and mesh
        // show_skeleton(viewer, CT, BET);
        // show_mesh_and_skeleton(viewer, VT, F, CT, BET);
        show_mesh(viewer, VT, F);

        frame = (frame + 1) % (transM.rows() / (3 * 20));
        anim_t += anim_t_dir;
        anim_t_dir *= (anim_t>=1.0 || anim_t<=0.0?-1.0:1.0);
    }

    return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods) {
  switch(key)
  {
    case ' ':
        viewer.core().is_animating = !viewer.core().is_animating;
        break;
    case 'S':       // show S elected handles
        show_handles(viewer, V, refH, handle_id);
        handle_id = (handle_id + 1) % (C.rows());
        break;
    case 'R':       // show R eference handles
        show_handles(viewer, V, refH, handle_id);
        handle_id = (handle_id + 1) % (C.rows());
        break;
    case 'W':       // show skinning W eight function
        show_skinning_weight_function(viewer, W, handle_id);
        handle_id = (handle_id + 1) % (C.rows());
        break;
    case 'V':       // show per-V ertex linear blending skinning
        // TODO show_per_vertex_linear_blending_skinning();
        break;
    case 'F':       // show per-F ace linear blending skinning
        // TODO show_per_vertex_linear_blending_skinning();
        break;
    case 'E':       // show example
        Eigen::MatrixXd EV;
        Eigen::MatrixXi EF;
        igl::readOFF("../data/GP/body.off", EV, EF);
        show_mesh(viewer, EV, EF);
        break;
    case 'U':       // show unposed example
        // TODO
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
	
    // load data
    load_mesh("../data/hand/hand.off", V, F);
    load_skeleton("../data/hand/hand.tgf", C, E);
    load_quat("../data/hand/hand-pose_quat.dmat", transQ);
    load_matrot("../data/hand/hand-pose_matrot.dmat", transM);
    load_handles("../data/hand/hand-handles.dmat", refH);

    // preprocess
    igl::directed_edge_parents(E, P);
    igl::column_to_quats(transQ, pose);
    save_handle_and_free_vertices(refH, ref_handle_vertices, ref_free_vertices);
    compute_root_vertices_and_positions(V, P, refH, ref_root_vertices, ref_non_root_vertices, ref_root_handle_positions);

    // show mesh and skeleton
    show_mesh_and_skeleton(viewer, V, F, C, E);

    // select handles
    // select_handles(C, E, V, H);
    // save_handle_and_free_vertices(H, handle_vertices, free_vertices);

    // compute per-vetex skinning weight function
    laplace_prefactor(V, F, Lw, lAff, lAfc, ref_handle_vertices, ref_free_vertices, laplace_solver);
    compute_per_vertex_skinning_weight_function(E, refH, lAfc, ref_handle_vertices, ref_free_vertices, laplace_solver, W);

    // compute per-face skinning weight function
    compute_per_face_skinning_weight_function(W, F, FW);
    compute_gradient_operator(V, F, G);
    compute_diagonal_weighting_matrix(V, F, D);
    RotationList fQ;
    poisson_prefactor(G, D, ref_root_vertices, ref_non_root_vertices, poisson_solver, pAff, pAfc);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    menu.callback_draw_viewer_menu = [&]() {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Visualization Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("visualize mesh", ImVec2(-1,0))) {
                show_mesh(viewer, V, F);
            }

            if (ImGui::Button("visualize skeleton", ImVec2(-1,0))) {
                show_skeleton(viewer, C, E);
            }

            if (ImGui::Button("visualize mesh & skeleton", ImVec2(-1,0))) {
                show_mesh_and_skeleton(viewer, V, F, C, E);
            }

            // if (ImGui::Button("visualize selected handles", ImVec2(-1,0))) {
            //     // TODO
            // }

            // if (ImGui::Button("visualize reference handles", ImVec2(-1,0))) {
            //     show_handles(viewer, V, refH, handle_id);
            // }
        }
    };

    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &key_down;
    viewer.core().is_animating = false;
    viewer.core().animation_max_fps = 30.;
    
    viewer.launch();
}