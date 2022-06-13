#include <iostream>
#include <math.h>
#include <igl/project_to_line.h>

using namespace std;

void compute_point_to_line_segment_distance (Eigen::RowVector3d v, Eigen::RowVector3d& s, Eigen::RowVector3d& d, double& sqrtd) {
        using namespace Eigen;

        RowVector3d p = v;
        RowVector3d sd = d - s;
        RowVector3d ds = s - d;
        RowVector3d sp = p - s;
        RowVector3d dp = p - d;

        double dotS = sd.dot(sp);
        double dotD = ds.dot(dp);
        if (dotS > 0 && dotD < 0) {
            // sqrtD[i] = dp.norm();
            sqrtd = 1000;
        } else if (dotS < 0 && dotD > 0) {
            // sqrtD[i] = sp.norm();
            sqrtd = 1000;
        } else {
            sqrtd = sd.cross(sp).norm() / sd.norm();
        }
    }

void compute_points_to_line_segment_distances (Eigen::MatrixXd& V, Eigen::RowVector3d& s, Eigen::RowVector3d& d,
    Eigen::VectorXd& sqrtD) {

    using namespace Eigen;

    sqrtD.resize(V.rows());
    
    for (int i = 0; i < V.rows(); i ++) {
        double sqrtd;
        compute_point_to_line_segment_distance (V.row(i), s, d, sqrtd);
        sqrtD[i] = sqrtd;
    }
}

// select handle vertices for joint [joint_id]
void select_handles_helper (const int joint_id, Eigen::MatrixXd& C, Eigen::MatrixXi& E, Eigen::MatrixXd& V, Eigen::VectorXi& H) {
    using namespace Eigen;

    RowVector3d s = C.row(E(joint_id, 0));
    RowVector3d d = C.row(E(joint_id, 1));

    // compute distance
    VectorXd sqrtD;
    compute_points_to_line_segment_distances(V, s, d, sqrtD);

    // find the distance to the nearest vertex
    double minDistance = sqrtD.minCoeff();

    // assign handle id to vertices
    double threshold = 2.5 * minDistance;
    for (int i = 0; i < V.rows(); i ++) {
        if (H[i] == -1 && sqrtD[i] < threshold) {
            H[i] = joint_id;
        } else if (H[i] > -1 && sqrtD[i] < threshold) {
            // TODO: how to deal with conflicting joint_id
            // RowVector3d ss = C.row(E(H[i], 0));
            // RowVector3d dd = C.row(E(H[i], 1));
            // double sqrtd;
            // compute_point_to_line_segment_distance(V.row(i), ss, dd, sqrtd);
            // if (sqrtd < sqrtD[i]) {
            //     H[i] = joint_id;
            // }
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
    Eigen::VectorXi& handle_vertices, Eigen::VectorXi& free_vertices) {
    
    int numFree = (H.array() == -1).cast<int>().sum();
    int numHandle = H.rows() - numFree;

    handle_vertices.resize(numHandle);
    free_vertices.resize(numFree);

    int count = 0;
    for (int i = 0; i < H.rows(); i ++) {
        if (H[i] > -1) {
            handle_vertices[count ++] = i;
        }
    }

    count = 0;
    for (int i = 0; i < H.rows(); i ++) {
        if (H[i] == -1) {
            free_vertices[count ++] = i;
        }
    }
}