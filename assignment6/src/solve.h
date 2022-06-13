#include <iostream>
#include <math.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

using namespace std;

Eigen::SparseMatrix<double> Lw, Aff, Afc;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> solver;

void prefactor(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& Lw,
    Eigen::SparseMatrix<double>& Aff, Eigen::SparseMatrix<double>& Afc,
    vector<int>& handle_vertices, vector<int>& free_vertices, 
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> solver) {

  igl::cotmatrix(V, F, Lw);

  igl::slice(Lw, free_vertices, free_vertices, Aff);
  igl::slice(Lw, free_vertices, handle_vertices, Afc);

  solver.compute(Aff);    // left hand side
}

void solve_laplacian(Eigen::SparseMatrix<double>& Afc,
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> solver,
    vector<int>& handle_vertices, const int k) {
    
    // compute right hand side
    

    return ;
}