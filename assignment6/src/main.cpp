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
#include "context_aware.h"
// #include "utils.h"

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// for task 8
Eigen::MatrixXd transM0, transM1, transM2, transM3;
Eigen::MatrixXd V0, V1, V2, V3;
Eigen::MatrixXd F0, F1, F2, F3;
vector<Eigen::MatrixXd> ET, UV;
Eigen::MatrixXd JC;

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
Eigen::MatrixXd W, refW, FW, fQG;
Eigen::SparseMatrix<double> Lw, G, D, lAff, lAfc, pAff, pAfc;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> laplace_solver;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> poisson_solver;

// for animation
RotationList pose;
int frame = 0;
int frame_dir = 1;
int numFrames;
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

// GUI
int animation_mode = 0;

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

void compute_absolute_transmat(Eigen::MatrixXd& transM, int frame, Eigen::MatrixXd& T) {
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

    return ;
}

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

// load and pre-compute for hand dataset
void load_hand () {
    // load data for task 3, 4, 5, 6, 7
    igl::readOFF("../data/hand/hand.off", V, F);
    load_skeleton("../data/hand/hand.tgf", C, E);
    load_quat("../data/hand/hand-pose_quat.dmat", transQ);
    load_matrot("../data/hand/hand-pose_matrot.dmat", transM);
    load_handles("../data/hand/hand-handles.dmat", refH);

    return ;
}

void hand_precompute () {
    numFrames = transM.rows() / (3 * (E.rows() - 1));
    igl::directed_edge_parents(E, P);
    igl::column_to_quats(transQ, pose);

    // compute per-vertex skinning weight function using reference handles
    save_handle_and_free_vertices(refH, ref_handle_vertices, ref_free_vertices);
    compute_root_vertices_and_positions(V, P, refH, ref_root_vertices, ref_non_root_vertices, ref_root_handle_positions);
    laplace_prefactor(V, F, Lw, lAff, lAfc, ref_handle_vertices, ref_free_vertices, laplace_solver);
    compute_per_vertex_skinning_weight_function(E, refH, lAfc, ref_handle_vertices, ref_free_vertices, laplace_solver, refW);

    // compute per-vertex skinning weight function using selected handles
    select_handles(C, E, V, H);
    save_handle_and_free_vertices(H, handle_vertices, free_vertices);
    compute_root_vertices_and_positions(V, P, H, root_vertices, non_root_vertices, root_handle_positions);
    laplace_prefactor(V, F, Lw, lAff, lAfc, handle_vertices, free_vertices, laplace_solver);
    compute_per_vertex_skinning_weight_function(E, H, lAfc, handle_vertices, free_vertices, laplace_solver, W);

    // compute per-face skinning weight function using reference handles
    compute_per_face_skinning_weight_function(refW, F, FW);
    compute_gradient_operator(V, F, G);
    compute_diagonal_weighting_matrix(V, F, D);
    poisson_prefactor(G, D, ref_root_vertices, ref_non_root_vertices, poisson_solver, pAff, pAfc);

    return ;
}

void load_context_aware () {
    // load data for task 8
    igl::readOBJ("../data/context-aware/reference.obj", V, F);
    load_skeleton("../data/context-aware/skeleton.tgf", C, E);
    load_matrot("../data/context-aware/all_frames.dmat", transM);
    load_handles("../data/context-aware/handles.dmat", refH);

    load_matrot("../data/context-aware/eg0.dmat", transM0);
    load_matrot("../data/context-aware/eg1.dmat", transM1);
    load_matrot("../data/context-aware/eg2.dmat", transM2);
    load_matrot("../data/context-aware/eg3.dmat", transM3);
    
    igl::readOBJ("../data/context-aware/eg0.obj", V0, F0);
    igl::readOBJ("../data/context-aware/eg1.obj", V1, F1);
    igl::readOBJ("../data/context-aware/eg2.obj", V2, F2);
    igl::readOBJ("../data/context-aware/eg3.obj", V3, F3);

    return ;
}

void context_aware_precompute () {
    numFrames = transM.rows() / (3 * (E.rows() - 1));
    igl::directed_edge_parents(E, P);
    igl::column_to_quats(transQ, pose);

    // compute per-vertex skinning weight function using reference handles
    save_handle_and_free_vertices(refH, ref_handle_vertices, ref_free_vertices);
    compute_root_vertices_and_positions(V, P, refH, ref_root_vertices, ref_non_root_vertices, ref_root_handle_positions);
    laplace_prefactor(V, F, Lw, lAff, lAfc, ref_handle_vertices, ref_free_vertices, laplace_solver);
    compute_per_vertex_skinning_weight_function(E, refH, lAfc, ref_handle_vertices, ref_free_vertices, laplace_solver, refW);

    // compute per-vertex skinning weight function using selected handles
    select_handles(C, E, V, H);
    save_handle_and_free_vertices(H, handle_vertices, free_vertices);
    compute_root_vertices_and_positions(V, P, H, root_vertices, non_root_vertices, root_handle_positions);
    laplace_prefactor(V, F, Lw, lAff, lAfc, handle_vertices, free_vertices, laplace_solver);
    compute_per_vertex_skinning_weight_function(E, H, lAfc, handle_vertices, free_vertices, laplace_solver, W);

    // unpose examples
    Eigen::MatrixXd ET0, ET1, ET2, ET3;
    compute_absolute_transmat(transM0, 0, ET0);
    compute_absolute_transmat(transM1, 0, ET1);
    compute_absolute_transmat(transM2, 0, ET2);
    compute_absolute_transmat(transM3, 0, ET3);
    ET.push_back(ET0);
    ET.push_back(ET1);
    ET.push_back(ET2);
    ET.push_back(ET3);
    compute_c(ET, JC);
    Eigen::MatrixXd UV0, UV1, UV2, UV3;
    unpose(V0, refW, ET0, UV0);
    unpose(V1, refW, ET1, UV1);
    unpose(V2, refW, ET2, UV2);
    unpose(V3, refW, ET3, UV3);
    UV.push_back(UV0);
    UV.push_back(UV1);
    UV.push_back(UV2);
    UV.push_back(UV3);
}

bool pre_draw (igl::opengl::glfw::Viewer & viewer) {
    using namespace Eigen;

    // compute absolute transformations
    MatrixXd T;
    RotationList vQ;
    vector<Vector3d> vT;
    MatrixXd CT;
    MatrixXi BET;
    MatrixXd VT;
    RotationList fQ;

    if(viewer.core().is_animating) {
        switch (animation_mode)
        {
            case -1:        // skeleton deformation
                // compute_absolute_transmat(transM, frame, T);
                compute_absolute_transquat(pose, anim_t, T);
                igl::deform_skeleton(C, E, T, CT, BET);
                show_skeleton(viewer, CT, BET);
                break;
            case 0:         // per-vertex linear blending skinning
                compute_absolute_transmat(transM, frame, T);
                // compute_absolute_transquat(pose, anim_t, T);
                compute_per_vertex_linear_blending_skinning(V, refW, T, VT);
                show_mesh(viewer, VT, F);
                break;
            case 1:         // dual quaternion skinning
                compute_absolute_RTquat(pose, anim_t, vQ, vT);
                compute_dual_quaternion_skinning(V, W, vQ, vT, VT);
                show_mesh(viewer, VT, F);
                break;
            case 2:         // per-face linear blending skinning
                compute_absolute_RTquat(pose, anim_t, vQ, vT);
                compute_per_face_transformation(FW, vQ, fQ);
                compute_face_deformation_gradient(fQ, fQG);
                compute_poisson_stitching(V, G, D, fQG, poisson_solver, pAfc, ref_root_vertices, ref_non_root_vertices, ref_root_handle_positions, VT);
                show_mesh(viewer, VT, F);
                break;
            case 3:         // per-vertex context-aware linear blending skinning
                compute_absolute_transmat(transM, frame, T);
                compute_per_vertex_context_aware_deformation (V, UV, refW, T, ET, JC, VT);
                show_mesh(viewer, VT, F);
                break;
        }
        // if (matrot) {
            // compute_absolute_transmat(transM, frame, T);
            // // cout << "frame: " << frame << endl;
        // } else {
            // compute_absolute_transquat(pose, anim_t, T);
            // compute_absolute_RTquat(pose, anim_t, vQ, vT);
        // }

        // skeleton deformation
        // MatrixXd CT;
        // MatrixXi BET;
        // igl::deform_skeleton(C, E, T, CT, BET);

        // per-vertex linear blending skinning
        // MatrixXd VT;
        // compute_per_vertex_linear_blending_skinning(V, W, T, VT);

        // dual quaternion skinning
        // compute_dual_quaternion_skinning(V, W, vQ, vT, VT);

        // per-face linear blending skinning
        // RotationList fQ;
        // compute_per_face_transformation(FW, vQ, fQ);
        // compute_face_deformation_gradient(fQ, fQG);
        // compute_poisson_stitching(V, G, D, fQG, poisson_solver, pAfc, ref_root_vertices, ref_non_root_vertices, ref_root_handle_positions, VT);

        // per-vertex context-aware deformation
        // compute_per_vertex_context_aware_deformation (V, UV, W, T, ET, JC, VT);

        // display deformed skeleton and mesh
        // show_skeleton(viewer, CT, BET);
        // show_mesh_and_skeleton(viewer, VT, F, CT, BET);
        // show_mesh(viewer, VT, F);

        if (frame == numFrames - 2) {
            frame_dir = -1;
        } else if (frame == 0) {
            frame_dir = 1;
        }
        frame += frame_dir;
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
    case 'H':       // load and precompute for hand dataset
        load_hand();
        hand_precompute();
        show_mesh_and_skeleton(viewer, V, F, C, E);
        break;
    case 'B':
        load_context_aware();       // load and precompute for context-aware (body) dataset
        context_aware_precompute();
        show_mesh_and_skeleton(viewer, V, F, C, E);
        break;
    case 'S':       // show S elected handles
        show_handles(viewer, V, H, handle_id);
        handle_id = (handle_id + 1) % (C.rows());
        break;
    case 'R':       // show R eference handles
        show_handles(viewer, V, refH, handle_id);
        handle_id = (handle_id + 1) % (C.rows());
        break;
    case 'W':       // show skinning W eight function computed on selected handles
        show_skinning_weight_function(viewer, W, handle_id);
        handle_id = (handle_id + 1) % (C.rows());
        break;
    case 'E':       // show skeleton deformation
        animation_mode = -1;
        break;
    case 'V':       // show per-V ertex linear blending skinning
        animation_mode = 0;
        break;
    case 'D':       // show dual quaternion skinning
        animation_mode = 1;
        break;
    case 'P':       // show per-face linear blending skinning
        animation_mode = 2;
        break;
    case 'Q':       // show per-vertex context-aware linear blending skinning
        animation_mode = 3;
        break;
    case 'U':       // show unposed example
        // MatrixXd& T, UV;
        // compute_absolute_transmat(transM, 0, T);
        // unpose(V1, refW, T, UV);
        // show_mesh(viewer, UV, F);
        break;
  }

  return true;
}

int main(int argc, char *argv[]) {
    Viewer viewer;
    
    std::string filename;
    if (argc == 2) {
        filename = std::string(argv[1]); // Mesh provided as command line argument
    }
    else {
        filename = std::string("../data/hand/hand.off"); // Default mesh
    }
    
    // load and show mesh
    if (filename.substr(filename.length() - 4) == ".off")
    {   
        igl::readOFF(filename, V, F);
        show_mesh(viewer, V, F);
    }
    else if (filename.substr(filename.length() - 4) == ".obj")
    {
        igl::readOBJ(filename, V, F);
        show_mesh(viewer, V, F);
    }
    else
    {
        std::cerr << "Extension unknown (must be '.off' or '.obj')\n";
        return false;
    }

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