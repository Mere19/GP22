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

// typedef Eigen::Triplet<double> Triplet;

void compute_T_sum (Eigen::MatrixXd& V, Eigen::MatrixXd& W, Eigen::MatrixXd& T, Eigen::MatrixXd& T_sum) {
    using namespace Eigen;

    int numVertices = V.rows();
    int numJoints = W.cols();

    T_sum.setZero(3 * numVertices, 4);
    for (int i = 0; i < numVertices; i ++) {
        for (int j = 0; j < numJoints; j ++) {
            T_sum.block(3 * i, 0, 3, 4) += W(i, j) * T.block(4 * j, 0, 4, 3).transpose();
        }
    }

    // cout << "[completed] compute T sum" << endl;

    return ;
}

// T: stack of 3 by 4 transformations
// UV: unposed vertices
void unpose (Eigen::MatrixXd V, Eigen::MatrixXd& W, Eigen::MatrixXd& T, Eigen::MatrixXd& UV) {
    using namespace Eigen;

    int numVertices = V.rows();

    MatrixXd T_sum;
    compute_T_sum(V, W, T, T_sum);

    UV.resize(numVertices, 3);
    for (int i = 0; i < numVertices; i ++) {
        MatrixXd R = T_sum.block(3 * i, 0, 3, 3);
        Vector3d t = T_sum.block(3 * i, 3, 3, 1);
        UV.row(i) = (R.inverse() * (V.row(i).transpose() - t)).transpose();
    }

    return ;
}

// void compute_per_vertex_displacements (Eigen::MatrixXd& V, Eigen::MatrixXd& UV, Eigen::MatrixXd& S) {
//     S = UV - V;

//     return ;
// }

// void compute_delta_vector (Eigen::MatrixXd& V, vector<Eigen::MatrixXd>& UV) {
//     // TODO
// }

double compute_rbf (double x) {
    return exp(- x * x / 2);
}

void compute_rbf_vector (Eigen::MatrixXd& T, vector<Eigen::MatrixXd>& ET, Eigen::VectorXd& vecRBF) {
    vecRBF.setZero(ET.size());
    for (int i = 0; i < ET.size(); i ++) {
        vecRBF[i] = compute_rbf((T - ET[i]).norm());
    }
}

// JC: JC.col(j) contains the J unknowns c(j, t) for example j
void compute_c (vector<Eigen::MatrixXd>& ET, Eigen::MatrixXd& JC) {
    using namespace Eigen;

    int numExamples = ET.size();

    JC.resize(numExamples, numExamples);

    for (int i = 0; i < numExamples; i ++) {
        // solve for JC.col(i)
        MatrixXd A(numExamples, numExamples);       // construct A
        for (int j = 0; j < numExamples; j ++) {
            if (j == i) {
                A.row(j).setOnes();
            } else {
                VectorXd vecRBF;
                compute_rbf_vector(ET[j], ET, vecRBF);
                A.row(j) = vecRBF;
            }
        }

        VectorXd b;     // construct b
        b.setZero(numExamples);
        b[i] = 1;

        VectorXd x = A.colPivHouseholderQr().solve(b);       // solve for x
        JC.col(i) = x;
    }
    
    return ;
}

void compute_a (Eigen::MatrixXd& T, vector<Eigen::MatrixXd>& ET, Eigen::MatrixXd& JC, Eigen::VectorXd& a) {
    using namespace Eigen;

    int numExamples = ET.size();

    a.resize(numExamples);

    VectorXd RBFvec;
    compute_rbf_vector(T, ET, RBFvec);
    for (int i = 0; i < numExamples; i ++) {
        a[i] = JC.col(i).dot(RBFvec);
    }
    a = a / a.sum();

    return ;
}

void compute_delta(Eigen::MatrixXd& V, vector<Eigen::MatrixXd>& UV, int i, Eigen::MatrixXd& delta) {
    using namespace Eigen;

    int numExamples = UV.size();

    delta.resize(numExamples, 3);
    for (int j = 0; j < numExamples; j ++) {
        delta.row(j) = UV[j].row(i) - V.row(i);
    }

    return ;
}

// TODO.
// V: rest pose mesh vertices
void compute_per_vertex_context_aware_deformation (Eigen::MatrixXd& V, vector<Eigen::MatrixXd>& UV, Eigen::MatrixXd& W,
Eigen::MatrixXd& T, vector<Eigen::MatrixXd>& ET, Eigen::MatrixXd& JC, Eigen::MatrixXd& VT) {
    using namespace Eigen;

    // standard LBS
    MatrixXd M;
    igl::lbs_matrix(V, W, M);
    VT = M * T;

    // add context-aware information
    // for (int i = 0; i < V.rows(); i ++) {
    //     VectorXd a;
    //     MatrixXd delta;
    //     compute_a(T, ET, JC, a);
    //     compute_delta(V, UV, i, delta);
    //     VT.row(i) = VT.row(i) + a.transpose() * delta;
    // }

    return ;
}