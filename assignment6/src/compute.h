#include <iostream>
#include <math.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/normalize_row_sums.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/lbs_matrix.h>
#include <igl/dqs.h>
#include <igl/grad.h>

#include "utils.h"

using namespace std;

void laplace_prefactor(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& Lw,
    Eigen::SparseMatrix<double>& Aff, Eigen::SparseMatrix<double>& Afc,
    Eigen::VectorXi& handle_vertices, Eigen::VectorXi& free_vertices,
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor>& solver) {

    igl::cotmatrix(V, F, Lw);

    igl::slice(Lw, free_vertices, free_vertices, Aff);
    igl::slice(Lw, free_vertices, handle_vertices, Afc);

    solver.compute(Aff);    // left hand side
}

void solve_laplace (Eigen::VectorXi& H, Eigen::SparseMatrix<double>& Afc, Eigen::VectorXi& handle_vertices, Eigen::VectorXi& free_vertices, const int k,
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor>& solver, Eigen::VectorXd& Wk) {
    using namespace Eigen;

    int numHandle = handle_vertices.size();
    int numFree = H.rows() - handle_vertices.size();

    Wk.resize(H.rows());
    VectorXd vc(numHandle);

    // compute right hand side -Afc * v
    for (int i = 0; i < numHandle; i ++) {
        int idx = handle_vertices[i];
        if (H[idx] == k) {
            vc[i] = 1;
        } else {
            vc[i] = 0;
        }
    }

    // solve for Aff
    Eigen::VectorXd vf = solver.solve(-Afc * vc);

    // concatenate and save result
    // int count = 0;
    // for (int i = 0; i < H.rows(); i ++) {
    //     if (H[i] == -1 && count < numFree) {
    //       Wk[i] = vf[count ++];
    //     }
    // }
    igl::slice_into(vf, free_vertices, 1, Wk);
    for (int i = 0; i < H.rows(); i ++) {
        if (H[i] == k) {
          Wk[i] = 1;
        } else if (H[i] > -1 && H[i] != k) {
          Wk[i] = 0;
        }
    }

    return ;
}

void compute_per_vertex_skinning_weight_function (Eigen::MatrixXi& E, Eigen::VectorXi& H, Eigen::SparseMatrix<double>& Afc,
    Eigen::VectorXi& handle_vertices, Eigen::VectorXi& free_vertices,
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor>& solver, Eigen::MatrixXd& W) {

    using namespace Eigen;

    int numHandle = handle_vertices.size();
    int numFree = H.rows() - handle_vertices.size();

    W.resize(H.rows(), E.rows());
    
    for (int k = 0; k < E.rows(); k ++) {
      Eigen::VectorXd Wk;
      solve_laplace(H, Afc, handle_vertices, free_vertices, k, solver, Wk);
      W.col(k) = Wk;
    }
}

void compute_per_face_skinning_weight_function (Eigen::MatrixXd& W, Eigen::MatrixXi& F, Eigen::MatrixXd& FW) {
    using namespace Eigen;

    FW.resize(F.rows(), W.cols());

    for (int i = 0; i < FW.rows(); i ++) {
        RowVector3i v = F.row(i);
        for (int j = 0; j < FW.cols(); j ++) {
            double wv0 = W(v[0], j);
            double wv1 = W(v[1], j);
            double wv2 = W(v[2], j);

            FW(i, j) = (wv0 + wv1 + wv2) / 3.0;
        }
    }

    return ;
}

void compute_per_face_transformation (Eigen::MatrixXd& FW, RotationList& vQ, RotationList& fQ) {
    using namespace Eigen;

    fQ.resize(FW.rows());   // initialize

    // compute the face rotation
    for (int f = 0; f < FW.rows(); f ++) {
        Vector3d sum_log;
        sum_log.setZero();

        for (int k = 0; k < FW.cols(); k ++) {
            AngleAxisd qk;
            qk = vQ[k].normalized();
            Vector3d log_qk = -qk.angle() * qk.axis();
            sum_log += FW(f, k) * log_qk;
        }

        AngleAxisd exp(sum_log.norm(), sum_log.normalized());
        Quaterniond qf;
        qf = exp;
        fQ[f] = qf;
    }

    return ;
}

void compute_gradient_operator(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& G) {
    using namespace Eigen;

    // compute
    igl::grad(V, F, G);

    // permute
    int fnum = F.rows();
    VectorXd rows_permuted, cols_permuted;
    rows_permuted.resize(3*fnum);
    cols_permuted.resize(V.rows());

    for (int i = 0, j = 0; i < 3*fnum, j < fnum; i +=3, j ++) {
      rows_permuted[i] = j;
      rows_permuted[i + 1] = j + fnum;
      rows_permuted[i + 2] = j + 2*fnum;
    }

    for (int i = 0; i < V.rows(); i ++) {
      cols_permuted[i] = i;
    }

    igl::slice(G, rows_permuted, cols_permuted, G);

    return ;
}

void compute_face_deformation_gradient (RotationList& fQ, Eigen::MatrixXd& fQG) {
    using namespace Eigen;

    fQG.resize(fQ.size() * 3, 3);

    for (int i = 0; i < fQ.size(); i ++) {
        Matrix3d grad;
        grad = fQ[i];
        fQG.block(i * 3, 0, 3, 3) = grad;
    }

    return ;
}

void compute_diagonal_weighting_matrix(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<double>& D) {
  using namespace Eigen;

  VectorXd dblA;

  igl::doublearea(V, F, dblA);
  VectorXd sglA = dblA / 2.0;
  VectorXd sglA_vec;
  sglA_vec.resize(3*F.rows());

  for (int i = 0, j = 0; i < 3*F.rows(), j < F.rows(); i += 3, j ++) {
    sglA_vec[i] = sglA[j];
    sglA_vec[i+1] = sglA[j];
    sglA_vec[i+2] = sglA[j];
  }

  igl::diag(sglA_vec, D);

  return;
}

void poisson_prefactor (Eigen::SparseMatrix<double>& G, Eigen::SparseMatrix<double>& D,
Eigen::VectorXi& handle_vertices, Eigen::VectorXi& free_vertices,
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor>& solver,
Eigen::SparseMatrix<double>& Aff, Eigen::SparseMatrix<double>& Afc) {
    Eigen::SparseMatrix<double> A = G.transpose() * D * G;   // compute G.T * D * G

    igl::slice(A, free_vertices, free_vertices, Aff);
    igl::slice(A, free_vertices, handle_vertices, Afc);

    solver.compute(Aff);    // left hand side
}

void compute_poisson_stitching (Eigen::MatrixXd& V, Eigen::SparseMatrix<double>& G, Eigen::SparseMatrix<double>& D, Eigen::MatrixXd& fQG,
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor>& solver, Eigen::SparseMatrix<double>& Afc,
Eigen::VectorXi& root_vertices, Eigen::VectorXi& non_root_vertices, Eigen::MatrixXd& root_handle_positions, Eigen::MatrixXd& VT) {
    using namespace Eigen;

    VT.resize(V.rows(), 3);   // initialize

    // compute free vertices positions
    MatrixXd M = G.transpose() * D * fQG;
    MatrixXd R1;
    Vector3d dim(0, 1, 2);
    igl::slice(M, non_root_vertices, 1, R1);
    MatrixXd R2 = Afc * root_handle_positions;
    MatrixXd fVT = solver.solve(R1 - R2);

    // concatenate and save the result
    // VectorXd root_vertices(root_handle_indices.data());
    // VectorXd non_root_vertices(non_root_handle_indices.data());
    igl::slice_into(root_handle_positions, root_vertices, 1, VT);
    igl::slice_into(fVT, non_root_vertices, 1, VT);
}

void compute_per_vertex_linear_blending_skinning (Eigen::MatrixXd& V, Eigen::MatrixXd& W, Eigen::MatrixXd& T, Eigen::MatrixXd& VT) {
    using namespace Eigen;

    // compute lbs matrix
    Eigen::MatrixXd M;
    igl::lbs_matrix(V, W, M);

    VT = M * T;

    return ;
}

void compute_dual_quaternion_skinning (Eigen::MatrixXd& V, Eigen::MatrixXd& W, RotationList& vQ, vector<Eigen::Vector3d>& vT, Eigen::MatrixXd& VT) {
    using namespace Eigen;

    VT.resize(V.rows(), 3);

    igl::dqs(V, W, vQ, vT, VT);

    return;
}