#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/jet.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <igl/sum.h>
#include <igl/speye.h>
#include <igl/bfs.h>
#include <igl/cotmatrix.h>
#include <igl/principal_curvature.h>
#include <imgui/imgui.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/doublearea.h>
#include <igl/cotmatrix.h>
#include <igl/adjacency_list.h>
#include <igl/massmatrix.h>
#include <igl/gaussian_curvature.h>
#include <igl/principal_curvature.h>
/*** insert any libigl headers here ***/

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #Vx3
Eigen::MatrixXd V;
// Face array, #Fx3
Eigen::MatrixXi F;
//Face normals #Fx3
Eigen::MatrixXd FN;
//Vertex normals #Vx3
Eigen::MatrixXd VN;

// Per-vertex uniform normal array, #Vx3
Eigen::MatrixXd N_uniform;
// Per-vertex area-weighted normal array, #Vx3
Eigen::MatrixXd N_area;
// Per-vertex mean-curvature normal array, #Vx3
Eigen::MatrixXd N_meanCurvature;
// Per-vertex PCA normal array, #Vx3
// K for K-ring
int K = 2;
Eigen::MatrixXd N_PCA;
// Per-vertex quadratic fitted normal array, #Vx3
Eigen::MatrixXd fromCanonical;
Eigen::MatrixXd toCanonical;
Eigen::MatrixXd N_quadraticFit;

// Per-vertex mean curvature, #Vx3
Eigen::VectorXd K_mean;
// Per-vertex Gaussian curvature, #Vx3
Eigen::VectorXd K_Gaussian;
// Per-vertex minimal principal curvature, #Vx3
Eigen::VectorXd K_min_principal;
// Per-vertex maximal principal curvature, #Vx3
Eigen::VectorXd K_max_principal;
// Per-vertex color array, #Vx3
Eigen::MatrixXd colors_per_vertex;

// number of iterations
int Niter = 20;
// lambda
double lamda = 0.001;
// time step
double dt = 1;
// Explicitely smoothed vertex array, #Vx3
Eigen::MatrixXd V_expLap;
// Implicitely smoothed vertex array, #Vx3
Eigen::MatrixXd V_impLap;
// Bilateral smoothed vertex array, #Vx3
Eigen::MatrixXd V_bilateral;


Eigen::SparseMatrix<double> uniform_laplacian(Eigen::MatrixXd V, Eigen::MatrixXi F) {
    // initialize laplacian
    Eigen::SparseMatrix<double> L(V.rows(), V.rows());
    L.setZero();

    // initialize list of non-zero elements (row, column, value)
    std::vector<Eigen::Triplet<double> > tripletList;

    // insert element to the list (row, column, value)
    for (int i = 0; i < F.rows(); i ++) {
        for (int j = 0; j < 3; j ++) {
            int j1 = (j+1)%3;
            int j2 = (j+2)%3;

            tripletList.push_back(Eigen::Triplet<double>(F(i, j1), F(i, j2), 1));
            tripletList.push_back(Eigen::Triplet<double>(F(i, j2), F(i, j1), 1));
            tripletList.push_back(Eigen::Triplet<double>(F(i, j1), F(i, j1), -1));
            tripletList.push_back(Eigen::Triplet<double>(F(i, j2), F(i, j2), -1));
        }
    }

    //construct matrix from the list
    L.setFromTriplets(tripletList.begin(), tripletList.end());

    return L;
}

std::vector<int> getAdaptiveNeighbours(int i, Eigen::MatrixXd V, std::vector<std::vector<int>> A, double sigma_c) {
    std::vector<int> neighbours;

    std::vector<bool> visited(V.rows(), false);
    std::queue<int> toBeVisited;
    visited[i] = true;
    toBeVisited.push(i);
    double radius = 2 * sigma_c;
    Eigen::RowVector3d vi = V.row(i);

    while (!toBeVisited.empty()) {
        int v = toBeVisited.front();
        toBeVisited.pop();
        neighbours.push_back(v);
        for (int j : A[v]) {
            if (visited[j] == false) {
                Eigen::RowVector3d vj = V.row(j);
                double dist = (vi - vj).norm();
                if (dist <= radius) {
                    toBeVisited.push(j);
                }
                visited[j] = true;
            }
        }
    }

    // cout << "num neighbours = " << neighbours.size() << endl;

    return neighbours;
}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing uniform vertex normals here:
        // store in N_uniform

        // initialize uniform norm
        N_uniform.resize(V.rows(), 3);

        // compute per face/vertex normal
        igl::per_face_normals(V, F, FN);
        igl::per_vertex_normals(V, F, VN);

        // compute vertex-face adjacency
        std::vector<std::vector<int>> VF;
        std::vector<std::vector<int>> VFi;
        igl::vertex_triangle_adjacency(V, F, VF, VFi);

        // compute uniform vertex normal for each vertex
        for (int i = 0; i < V.rows(); i ++) {
            Eigen::RowVector3d vertex_normal(0.0, 0.0, 0.0);
            std::vector<int> faces = VF[i];
            for (int f : faces) {
                vertex_normal += FN.row(f);
            }
            vertex_normal /= faces.size();

            // Use igl::per_vertex_normals to orient your normals consistently.
            if (vertex_normal * VN.row(i).transpose() < 0) {
                N_uniform.row(i) = -vertex_normal;
            } else {
                N_uniform.row(i) = vertex_normal;
            }
        }

        // Set the viewer normals.
        N_uniform.rowwise().normalize();
        viewer.data().set_normals(N_uniform);
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing area-weighted vertex normals here:
        // store in N_area

        // initialize area norm
        N_area.resize(V.rows(), 3);

        // compute per face/vertex normal
        igl::per_face_normals(V, F, FN);
        igl::per_vertex_normals(V, F, VN);

        // compute vertex-face adjacency
        std::vector<std::vector<int>> VF;
        std::vector<std::vector<int>> VFi;
        igl::vertex_triangle_adjacency(V, F, VF, VFi);

        // compute double face areas
        Eigen::VectorXd doubA;
        igl::doublearea(V, F, doubA);

        // compute uniform vertex normal for each vertex
        for (int i = 0; i < V.rows(); i ++) {
            Eigen::RowVector3d vertex_normal(0.0, 0.0, 0.0);
            std::vector<int> faces = VF[i];
            double total_weight = 0;
            for (int f : faces) {
                total_weight += doubA[f];
                vertex_normal += doubA[f] * FN.row(f);
            }
            vertex_normal /= total_weight;

            // Use igl::per_vertex_normals to orient your normals consistently.
            if (vertex_normal * VN.row(i).transpose() < 0) {
                N_area.row(i) = -vertex_normal;
            } else {
                N_area.row(i) = vertex_normal;
            }
        }

        // Set the viewer normals.
        N_area.rowwise().normalize();
        viewer.data().set_normals(N_area);
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing mean-curvature vertex normals here:
        // store in N_meanCurvature

        // initialize mean curvature norm
        N_meanCurvature.resize(V.rows(), 3);

        // compute cotangent matrix (Laplacian not scaled by area)
        Eigen::SparseMatrix<double> L;
        igl::cotmatrix(V, F, L);

        // compute cotangent-weighted discrete laplacian operator
        N_meanCurvature = L * V;

        // Use igl::per_vertex_normals to orient your normals consistently.
        for (int i = 0; i < V.rows(); i ++) {
            if (N_meanCurvature.row(i) * VN.row(i).transpose() < 0) {
                N_meanCurvature.row(i) = -N_meanCurvature.row(i);
            }
        }
    
        // Set the viewer normals.
        N_meanCurvature.rowwise().normalize();
        viewer.data().set_normals(N_meanCurvature);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing PCA vertex normals here:
        // store in N_PCA

        // initialize N_PCA
        N_PCA.resize(V.rows(), 3);

        // initialize and compute adjacency lists for each vertex
        std::vector<std::vector<int>> A;
        igl::adjacency_list(F, A);

        // compute PCA normal for each vertex
        for (int i = 0; i < V.rows(); i ++) {
            // BFS: initialize and compute k-ring neighbours for each vertex
            std::queue<int> to_be_visited;
            std::vector<int> visited(V.rows(), 0);
            std::vector<int> k_neighbours;

            int idx;
            visited[i] = 1;
            to_be_visited.push(i);
            k_neighbours.push_back(i);

            for (int j = 0; j < K; j ++) {
                std::queue<int> next_to_be_visited;
                while (!to_be_visited.empty()) {
                    idx = to_be_visited.front();
                    to_be_visited.pop();
                    for (int v : A[idx]) {
                        if (visited[v] == 0) {
                            visited[v] = 1;
                            next_to_be_visited.push(v);
                            k_neighbours.push_back(v);
                        }
                    }
                }
                to_be_visited = next_to_be_visited;
            }

            // compute k-ring neighbours' coordinates
            Eigen::MatrixXd X(k_neighbours.size(), 3);
            for (int j = 0; j < k_neighbours.size(); j ++) {
                int v = k_neighbours[j];
                X.row(j) = V.row(v);
            }

            // compute Y and S = Y * Y.t
            Eigen::MatrixXd Y = (X.rowwise() - X.colwise().mean()).transpose();
            Eigen::MatrixXd S = Y * Y.transpose();

            // compute eigenvector with the smallest eigenvalue for S
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
            Eigen::MatrixXd eigenvectors = es.eigenvectors();
            N_PCA.row(i) = eigenvectors.col(0).transpose();
        }

        // Use igl::per_vertex_normals to orient your normals consistently.
        for (int i = 0; i < V.rows(); i ++) {
            if (N_PCA.row(i) * VN.row(i).transpose() < 0) {
                N_PCA.row(i) = -N_PCA.row(i);
            }
        }

        // Set the viewer normals.
        N_PCA.rowwise().normalize();
        viewer.data().set_normals(N_PCA);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing quadratic fitted vertex normals here:
        // store in N_quadraticFit

        // initialize N_quadraticFit
        N_quadraticFit.resize(V.rows(), 3);

        // initialize and compute adjacency lists for each vertex
        std::vector<std::vector<int>> A;
        igl::adjacency_list(F, A);

        // compute quadratic fitting normal for each vertex
        for (int i = 0; i < V.rows(); i ++) {
            // BFS: initialize and compute k-ring neighbours for each vertex
            std::queue<int> to_be_visited;
            std::vector<int> visited(V.rows(), 0);
            std::vector<int> k_neighbours;

            int idx;
            visited[i] = 1;
            to_be_visited.push(i);
            k_neighbours.push_back(i);

            for (int j = 0; j < K; j ++) {
                std::queue<int> next_to_be_visited;
                while (!to_be_visited.empty()) {
                    idx = to_be_visited.front();
                    to_be_visited.pop();
                    for (int v : A[idx]) {
                        if (visited[v] == 0) {
                            visited[v] = 1;
                            next_to_be_visited.push(v);
                            k_neighbours.push_back(v);
                        }
                    }
                }
                to_be_visited = next_to_be_visited;
            }

            // compute k-ring neighbours' coordinates
            Eigen::MatrixXd X(k_neighbours.size(), 3);
            for (int j = 0; j < k_neighbours.size(); j ++) {
                int v = k_neighbours[j];
                X.row(j) = V.row(v);
            }
            
            // compute Y and S = Y * Y.t
            Eigen::MatrixXd X_norm = X.rowwise() - X.colwise().mean();
            Eigen::MatrixXd Y = X_norm.transpose();
            Eigen::MatrixXd S = Y * Y.transpose();

            // compute eigenvector with the smallest eigenvalue for S
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
            Eigen::MatrixXd eigenvectors = es.eigenvectors();
            double detS = eigenvectors.determinant();
            
            // compute change of basis matrix P (B' to B)
            if (detS < 0) {
                toCanonical = -eigenvectors;
            } else {
                toCanonical = eigenvectors;
            }
            // compute change of basis matrix P_invserse
            fromCanonical = toCanonical.inverse();

            // compute k-ring neighbours' coordinates in basis B'
            Eigen::MatrixXd X_norm_B_ = X_norm * fromCanonical.transpose();

            // compute W (quadratic bi-variate polynomial, u^2, v^2, uv, u, v, 1)
            Eigen::MatrixXd W(X_norm_B_.rows(), 6);
            for (int j = 0; j < X_norm_B_.rows(); j ++) {
                Eigen::RowVector3d x = X_norm_B_.row(j);
                // u^2
                W(j, 0) = x[1] * x[1];
                // v^2
                W(j, 1) = x[2] * x[2];
                // uv
                W(j, 2) = x[1] * x[2];
                // u
                W(j, 3) = x[1];
                // v
                W(j, 4) = x[2];
                // 1
                W(j, 5) = 1;
            }

            // compute h
            Eigen::VectorXd h = X_norm_B_.col(0);

            // solve Wb = h
            Eigen::VectorXd b = (W).colPivHouseholderQr().solve(h);

            // compute tangent vectors
            Eigen::RowVector3d vb_ = (V.row(i) - X.colwise().mean()) * fromCanonical.transpose();
            Eigen::RowVector3d t1(2*b[0]*vb_[1] + b[2] * vb_[2] + b[3], 1, 0);
            Eigen::RowVector3d t2(2*b[1]*vb_[2] + b[2] * vb_[1] + b[4], 0, 1);

            // compute normal
            N_quadraticFit.row(i) = t1.cross(t2) * toCanonical.transpose();
        }

        // Use igl::per_vertex_normals to orient your normals consistently.
        for (int i = 0; i < V.rows(); i ++) {
            if (N_quadraticFit.row(i) * VN.row(i).transpose() < 0) {
                N_quadraticFit.row(i) = -N_quadraticFit.row(i);
            }
        }

        // Set the viewer normals.
        N_quadraticFit.rowwise().normalize();
        viewer.data().set_normals(N_quadraticFit);
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete mean curvature:
        // store in K_mean

        Eigen::MatrixXd PD1, PD2;
        igl::principal_curvature(V, F, PD1, PD2, K_min_principal, K_max_principal);

        // compute mean curvature
        K_mean = 0.5 * (K_min_principal + K_max_principal);

        // For visualization, better to normalize the range of K_mean with the maximal and minimal curvatures.
        // store colors in colors_per_vertex
        Eigen::VectorXd K_mean_norm = K_mean / (K_mean.maxCoeff() - K_mean.minCoeff());
        igl::jet(K_mean_norm, true, colors_per_vertex);

        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '7') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete Gaussian curvature:
        // store in K_Gaussian

        // compute gaussian curvature
        igl::gaussian_curvature(V, F, K_Gaussian);

        // Compute mass matrix
        Eigen::SparseMatrix<double> M,Minv;
        igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
        igl::invert_diag(M,Minv);

        // Divide by area to get integral average
        K_Gaussian = (Minv*K_Gaussian).eval();

        // For visualization, better to normalize the range of K_Gaussian with the maximal and minimal curvatures.
        // store colors in colors_per_vertex
        Eigen::VectorXd K_Gaussian_norm = K_Gaussian / (K_Gaussian.maxCoeff() - K_Gaussian.minCoeff());
        igl::jet(K_Gaussian_norm, true, colors_per_vertex);

        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '8') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete minimal principal curvature:
        // store in K_min_principal
        
        Eigen::MatrixXd PD1, PD2;
        igl::principal_curvature(V, F, PD1, PD2, K_min_principal, K_max_principal);
       
        // For visualization, better to normalize the range of K_min_principal with the maximal and minimal curvatures.
        // store colors in colors_per_vertex
        Eigen::VectorXd K_min_principal_norm = K_min_principal / (K_min_principal.maxCoeff() - K_min_principal.minCoeff());
        igl::jet(K_min_principal_norm, true, colors_per_vertex);

        // Uncomment the code below to draw a blue segment parallel to the minimal curvature direction, 
        
        // const double avg = igl::avg_edge_length(V,F);
        // Eigen::Vector3d blue(0.2,0.2,0.8);
        // viewer.data().add_edges(V + PD_min*avg, V - PD_min*avg, blue);
        
        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '9') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete maximal principal curvature:
        // store in K_max_principal

        Eigen::MatrixXd PD1, PD2;
        igl::principal_curvature(V, F, PD1, PD2, K_min_principal, K_max_principal);

        // For visualization, better to normalize the range of K_max_principal with the maximal and minimal curvatures
        // store colors in colors_per_vertex
        Eigen::VectorXd K_max_principal_norm = K_max_principal / (K_max_principal.maxCoeff() - K_max_principal.minCoeff());
        igl::jet(K_max_principal_norm, true, colors_per_vertex);
        
        // Uncomment the code below to draw a red segment parallel to the maximal curvature direction
        
        // const double avg = igl::avg_edge_length(V,F);
        // Eigen::Vector3d red(0.8,0.2,0.2);
        // viewer.data().add_edges(V + PD_max*avg, V - PD_max*avg, red);
        
        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == 'E') {
        bool cotangent = false;
        // Add your code for computing explicit Laplacian smoothing here:
        // store the smoothed vertices in V_expLap
        V_expLap = V;

        // compute identity matrix
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(V.rows(), V.rows());

        // compute uniform-weighted Laplacian
        Eigen::SparseMatrix<double> L_uniform = uniform_laplacian(V, F);

        // compute cotangent-weighted Laplacian
        Eigen::SparseMatrix<double> L_cotangent;
        igl::cotmatrix(V, F, L_cotangent);

        // initialize mass matrix
        Eigen::SparseMatrix<double> M, Minv;

        Eigen::SparseMatrix<double> L;
        if (cotangent) {
            L = L_cotangent;
        } else {
            L = L_uniform;
        }

        // iteratively run the smoothing process
        for (int i = 0; i < Niter; i ++) {
            if (cotangent) {
                // compute mass matrix
                igl::massmatrix(V_expLap, F, igl::MASSMATRIX_TYPE_VORONOI, M);
                igl::invert_diag(M, Minv);
                // update
                V_expLap += dt * lamda * Minv * L * V_expLap;
            } else {
                // update
                V_expLap += dt * lamda * L * V_expLap;
            }
        }

        // Set the smoothed mesh
        viewer.data().clear();
        viewer.data().set_mesh(V_expLap, F);
    }

    if (key == 'D'){
        bool cotangent = true;
        // Implicit smoothing for comparison
        // store the smoothed vertices in V_impLap
        V_impLap = V;

        // compute identity matrix
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(V.rows(), V.rows());

        // compute uniform-weighted Laplacian
        Eigen::SparseMatrix<double> L_uniform = uniform_laplacian(V, F);

        // compute cotangent-weighted Laplacian
        Eigen::SparseMatrix<double> L_cotangent;
        igl::cotmatrix(V, F, L_cotangent);

        // initialize mass matrix
        Eigen::SparseMatrix<double> M, Minv;

        Eigen::SparseMatrix<double> L;
        if (cotangent) {
            L = L_cotangent;
        } else {
            L = L_uniform;
        }

        // iteratively run the smoothing process
        for (int i = 0; i < Niter; i ++) {
            if (cotangent) {
                // compute mass matrix
                igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
                // Solve (M-lamda*L) V = M*V
                const auto & S = (M - lamda*L);
                Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
                assert(solver.info() == Eigen::Success);
                V_impLap = solver.solve(M*V_impLap).eval();

                // Compute centroid and subtract (also important for numerics)
                Eigen::VectorXd dblA;
                igl::doublearea(V_impLap, F,dblA);
                double area = 0.5*dblA.sum();
                Eigen::MatrixXd BC;
                igl::barycenter(V_impLap, F, BC);
                Eigen::RowVector3d centroid(0,0,0);
                for(int i = 0;i<BC.rows();i++) {
                    centroid += 0.5 * dblA(i) / area * BC.row(i);
                }
                V_impLap.rowwise() -= centroid;
            } else {
                // update
            }
        }

        // Set the smoothed mesh
        viewer.data().clear();
        viewer.data().set_mesh(V_impLap, F);
    }

    if (key == 'B') {
        // Add your code for computing bilateral smoothing here:
        // store the smoothed vertices in V_bilateral
        // be care of the sign mistake in the paper
        // use v' = v - n * (sum / normalizer) to update
        V_bilateral = V;

        // other initialization

        for (int k = 0; k < 5; k ++) {
            cout << "k = " << k << endl;
            // initialize and compute adjacency lists for each vertex
            std::vector<std::vector<int>> A;
            igl::adjacency_list(F, A);

            // recompute normals
            igl::per_face_normals(V_bilateral, F, FN);
            igl::per_vertex_normals(V_bilateral, F, VN);

            // initialize area norm
            N_area.resize(V.rows(), 3);

            // compute per face/vertex normal
            igl::per_face_normals(V, F, FN);
            igl::per_vertex_normals(V, F, VN);

            // compute vertex-face adjacency
            std::vector<std::vector<int>> VF;
            std::vector<std::vector<int>> VFi;
            igl::vertex_triangle_adjacency(V, F, VF, VFi);

            // compute double face areas
            Eigen::VectorXd doubA;
            igl::doublearea(V, F, doubA);

            // compute uniform vertex normal for each vertex
            for (int i = 0; i < V.rows(); i ++) {
                Eigen::RowVector3d vertex_normal(0.0, 0.0, 0.0);
                std::vector<int> faces = VF[i];
                double total_weight = 0;
                for (int f : faces) {
                    total_weight += doubA[f];
                    vertex_normal += doubA[f] * FN.row(f);
                }
                vertex_normal /= total_weight;

                // Use igl::per_vertex_normals to orient your normals consistently.
                if (vertex_normal * VN.row(i).transpose() < 0) {
                    N_area.row(i) = -vertex_normal;
                } else {
                    N_area.row(i) = vertex_normal;
                }
            }
            N_area.rowwise().normalize();

            for (int i = 0; i < V.rows(); i ++) {
                Eigen::RowVector3d vi = V_bilateral.row(i);
                // determine sigma_c
                double sigma_c = 10;
                for (int j : A[i]) {
                    Eigen::RowVector3d vj = V_bilateral.row(j);
                    double dist = (vi - vj).norm();
                    if (dist < sigma_c) {
                        sigma_c = dist;
                    }
                }
                // get adaptive neighbours
                std::vector<int> neighbours = getAdaptiveNeighbours(i, V_bilateral, A, sigma_c);

                // parameter tuning
                double t, h, wc, ws, sum = 0, normalizer = 0;
                double sigma_s_sqrd = 0, h_mean = 0;
                for (int j : neighbours) {
                    Eigen::RowVector3d v_qi = vi - V_bilateral.row(j);
                    t = v_qi.norm();
                    h = v_qi * N_area.row(i).transpose();
                    h_mean += h;
                }
                h_mean /= neighbours.size();

                for (int j : neighbours) {
                    Eigen::RowVector3d v_qi = vi - V_bilateral.row(j);
                    h = v_qi * N_area.row(i).transpose();
                    sigma_s_sqrd += (h - h_mean) * (h - h_mean);
                }
                sigma_s_sqrd /= A[i].size();

                for (int j : neighbours) {
                    Eigen::RowVector3d v_qi = vi - V_bilateral.row(j);
                    t = v_qi.norm();
                    h = v_qi * N_area.row(i).transpose();
                    wc = exp(-t*t/(2*sigma_c*sigma_c));
                    ws = exp(-h*h/(2*sigma_s_sqrd));
                    sum += (wc * ws) * h;
                    normalizer += wc * ws;
                }

                V_bilateral.row(i) = vi - (sum / normalizer) * N_area.row(i);
            }
        }

        // Set the smoothed mesh
        viewer.data().clear();
        viewer.data().set_mesh(V_bilateral, F);
    }


    return true;
}

bool load_mesh(Viewer& viewer,string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    if (filename.substr(filename.length() - 4) == ".off")
    {
        igl::readOFF(filename, V, F);
    }
    else if (filename.substr(filename.length() - 4) == ".obj")
    {
        igl::readOBJ(filename, V, F);
    }
    else
    {
        std::cerr << "Extension unknown (must be '.off' or '.obj')\n";
        return false;
    }
    viewer.data().clear();
    viewer.data().set_mesh(V,F);
    viewer.data().compute_normals();
    viewer.core().align_camera_center(V, F);
    return true;
}

int main(int argc, char *argv[]) {
    // Show the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    
    std::string filename;
    if (argc == 2) {
        filename = std::string(argv[1]);
    }
    else {
        filename = std::string("../data/bumpy-cube.obj");
    }
    load_mesh(viewer,filename,V,F);

    callback_key_down(viewer, '1', 0);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
    viewer.launch();
}
