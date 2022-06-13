#include <iostream>
#include <math.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include <igl/colormap.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

int mesh_id, skeleton_id;

// show mesh in the viewer
void show_mesh (Viewer& viewer, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    mesh_id = viewer.append_mesh();
    viewer.selected_data_index = 0;
    viewer.data(mesh_id).set_mesh(V, F);
    viewer.data().compute_normals();
    viewer.core().align_camera_center(V, F);

    return ;
}

// show skeleton in the viewer
void show_skeleton (Viewer& viewer, Eigen::MatrixXd& C, Eigen::MatrixXi& E) {
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
    // viewer.data().clear();
    skeleton_id = viewer.append_mesh();
    viewer.data(skeleton_id).set_points(C, ptColors);
    viewer.data(skeleton_id).set_mesh(cyndV, cyndF);
    viewer.data().set_colors(cyndC);
    viewer.data().compute_normals();

    return ;
}

// 1, write down what I am doing
// 2, write down every possible subtask that I can come up with
// 3, cross the ones that do not work / finished

void show_mesh_and_skeleton (Viewer& viewer,
    Eigen::MatrixXd& V, Eigen::MatrixXi& F,
    Eigen::MatrixXd& C, Eigen::MatrixXi& E) {

    viewer.data().clear();

    show_mesh(viewer, V, F);        // show mesh
    show_skeleton(viewer, C, E);        // show skeleton
}

void show_handles (Viewer& viewer, Eigen::MatrixXd& V, Eigen::VectorXi& H, const int handle_id) {
    Eigen::MatrixXd handleColor;
    handleColor.resize(V.rows(), 3);

    Eigen::VectorXi selectedH(V.rows());
    for (int i = 0; i < V.rows(); i ++) {
        if (H[i] == handle_id) {
            selectedH[i] = handle_id;
        } else {
            selectedH[i] = -1;
        }
    }

    igl::colormap(igl::COLOR_MAP_TYPE_PLASMA, selectedH, true, handleColor);

    viewer.data(mesh_id).set_colors(handleColor);

    return ;
}

// W: #V by K, skinning weight function
// k: handle id
void show_skinning_weight_function (Viewer& viewer, Eigen::VectorXd& W, const int k) {
    using namespace Eigen;

    VectorXd ones;

    Eigen::MatrixXd weightColor;
    weightColor.resize(W.rows(), 3);
    ones.setOnes(W.rows());
    weightColor.col(0) = W;
    weightColor.col(1) = W;
    weightColor.col(2) = ones - W;
    viewer.data(mesh_id).set_colors(weightColor);
}