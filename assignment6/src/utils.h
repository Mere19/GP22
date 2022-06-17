#include <iostream>
#include <math.h>

typedef
  std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>>
  RotationList;

using namespace std;

// forward kinematics for rotation matrices dQ
void forward_kinematics(
    const Eigen::MatrixXd & C,
    const Eigen::MatrixXi & E,
    const Eigen::VectorXi & P,
    const Eigen::MatrixXd & dQ,
    std::vector<Eigen::Matrix3d> & vQ,
    std::vector<Eigen::Vector3d> & vT) {
    
    using namespace std;
    using namespace Eigen;

    const int dim = C.cols();
    const int m = E.rows();
    assert(m == P.rows());
    vector<bool> computed(m,false);
    vQ.resize(m);
    vT.resize(m);

    // Dynamic programming
    function<void (int) > fk_helper = [&] (int b) {
        if(!computed[b]) {
            if(P(b) < 0) {
                // base case for roots
                vQ[b] = dQ.block(b * dim, 0, dim, dim);
                const Vector3d r = C.row(E(b, 0)).transpose();
                vT[b] = r - vQ[b] * r;
            } else {
                // otherwise first compute parent's
                const int p = P(b);
                fk_helper(p);
                MatrixXd m = dQ.block(b * dim, 0, dim, dim);
                vQ[b] = vQ[p] * m;
                const Vector3d r = C.row(E(b, 0)).transpose();
                vT[b] = vT[p] - vQ[b] * r + vQ[p] * r;
            }
            computed[b] = true;
        }
    };

    for(int b = 0;b<m;b++) {
        fk_helper(b);
    }
}