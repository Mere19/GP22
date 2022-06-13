#include <iostream>
#include <math.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/normalize_row_sums.h>
#include <igl/slice_into.h>

using namespace std;

void prefactor(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& Lw,
    Eigen::SparseMatrix<double>& Aff, Eigen::SparseMatrix<double>& Afc,
    Eigen::VectorXi& handle_vertices, Eigen::VectorXi& free_vertices,
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor>& solver) {
      cout << "prefactoring ..." << endl;
  igl::cotmatrix(V, F, Lw);

  igl::slice(Lw, free_vertices, free_vertices, Aff);
  igl::slice(Lw, free_vertices, handle_vertices, Afc);

  solver.compute(Aff);    // left hand side
}

void solve_laplace (Eigen::VectorXi& H, Eigen::SparseMatrix<double>& Afc, Eigen::VectorXi& handle_vertices, const int k,
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor>& solver, Eigen::VectorXd& W) {
    using namespace Eigen;

    int numHandle = handle_vertices.size();
    int numFree = H.rows() - handle_vertices.size();

    W.resize(H.rows());
    VectorXd vc(handle_vertices.size());

    // compute right hand side -Afc * v
    for (int i = 0; i < numHandle; i ++) {
        int idx = handle_vertices[i];
        if (H[idx] == k) {
            vc[i] = 1;
        } else {
            vc[i] = 0;
        }
    }

    Eigen::VectorXd vf = solver.solve(-Afc * vc);

    // TODO: concatenate and save result
    int count = 0;
    for (int i = 0; i < H.rows(); i ++) {
        if (H[i] == -1 && count < numFree) {
          W[i] = vf[count ++];
        }
    }
    for (int i = 0; i < H.rows(); i ++) {
        if (H[i] == k) {
          W[i] = 1;
        } else if (H[i] > -1 && H[i] != k) {
          W[i] = 0;
        }
    }

    // TODO: sanity check
    // if (igl::normalize_row_sums())

    return ;
}

// void compute_skinning_weight_function (Eigen::VectorXi& H, Eigen::SparseMatrix<double>& Afc, Eigen::VectorXi& handle_vertices,
//     Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor>& solver, Eigen::MatrixXd& W) {

//     W.resize(H.rows(), handle_vertices.size());
    
//     for (int i = 0; i < H.rows(); i ++) {

//     }
// }