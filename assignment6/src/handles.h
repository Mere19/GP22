#include <iostream>
#include <math.h>
#include <igl/project_to_line.h>

using namespace std;

void compute_points_to_line_segment_distances(Eigen::MatrixXd& V, Eigen::RowVector3d& s, Eigen::RowVector3d& d,
    Eigen::VectorXd& sqrtD) {

    using namespace Eigen;

    sqrtD.resize(V.rows());
    
    for (int i = 0; i < V.rows(); i ++) {
        RowVector3d p = V.row(i);
        RowVector3d sd = d - s;
        RowVector3d ds = s - d;
        RowVector3d sp = p - s;
        RowVector3d dp = p - d;
        double dotS = sd.dot(sp);
        double dotD = ds.dot(dp);
        if (dotS > 0 && dotD < 0) {
            // sqrtD[i] = dp.norm();
            sqrtD[i] = 1000;
        } else if (dotS < 0 && dotD > 0) {
            // sqrtD[i] = sp.norm();
            sqrtD[i] = 1000;
        } else {
            double temp = dotD / ds.norm();
            sqrtD[i] = sd.cross(sp).norm() / sd.norm();
        }
    }
}

// select handle vertices for joint [joint_id]
void select_handles_helper (const int joint_id, Eigen::MatrixXd& C, Eigen::MatrixXi& E, Eigen::MatrixXd& V, Eigen::VectorXi& H) {
    using namespace Eigen;

    RowVector3d s = C.row(E(joint_id, 1));
    RowVector3d d = C.row(E(joint_id, 0));

    // compute distance
    VectorXd sqrtD;
    compute_points_to_line_segment_distances(V, s, d, sqrtD);

    // find the distance to the nearest vertex
    double minDistance = sqrtD.minCoeff();

    // assign handle id to vertices
    double threshold = 3 * minDistance;
    for (int i = 0; i < V.rows(); i ++) {
        if (H[i] == -1 && sqrtD[i] < threshold) {
            H[i] = joint_id;
        } else if (H[i] ) {
            // TODO
        }
    }
}

// select for each joint (bone edge) the handle vertices
void select_handles(Eigen::MatrixXd& C, Eigen::MatrixXi& E, Eigen::MatrixXd& V, Eigen::VectorXi& H) {
    H.resize(V.rows());
    for (int i = 0; i < H.rows(); i ++) {
        H[i] = -1;
    }

    for (int i = 0; i < E.rows(); i ++) {
        select_handles_helper(i, C, E, V, H);
    }
}

// save indices for handle vertices and free vertices
void save_handle_and_free_vertices (Eigen::VectorXi& H,
    vector<int> handle_vertices, vector<int> free_vertices) {

    for (int i = 0; i < H.rows(); i ++) {
        if (H[i] == -1) {
            free_vertices.push_back(i);
        } else {
            handle_vertices.push_back(i);
        }
    }
}