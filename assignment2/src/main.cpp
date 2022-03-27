#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/bounding_box_diagonal.h>

#include <math.h>
#include <sys/time.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// for timing
struct timeval startTime;
struct timeval endTime;

// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Normals evaluated via PCA method, #P x3
Eigen::MatrixXd NP;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
int polyDegree = 1;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 80;

// Parameter: grid resolution (resolution - 1 is multiple of 4)
int resolution = 21;

// Intermediate result: grid points, at which the imlicit function will be evaluated, #G x3
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

// create_spatial_index()
int resolutionSpatialIndex;
double dimUnitX, dimUnitY, dimUnitZ;
int minX, minY, minZ, maxX, maxY, maxZ;
std::vector<std::vector<std::vector<std::vector<int>>>> gridToVertices;

// epsilon
double diagonal, global_epsilon;

// Functions
void createGrid();
void evaluateImplicitFunc();
void evaluateImplicitFunc_PolygonSoup();
void getLines();
void pcaNormal();
bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers);

// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid()
{
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines.resize(0, 6);
    grid_values.resize(0);
    V.resize(0, 3);
    F.resize(0, 3);
    FN.resize(0, 3);

    // Grid bounds: axis-aligned bounding box
    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();

    // Bounding box dimensions
    Eigen::RowVector3d dim = (bb_max - bb_min) * 1.25;
    Eigen::RowVector3d extraSpace = dim - (bb_max - bb_min);

    // Grid spacing
    const double dx = (dim[0]) / (double)(resolution - 1);
    const double dy = (dim[1]) / (double)(resolution - 1);
    const double dz = (dim[2]) / (double)(resolution - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolution * resolution * resolution, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min - (extraSpace / 2) + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }
}

// Function for explicitly evaluating the implicit function for a sphere of
// radius r centered at c : f(p) = ||p-c|| - r, where p = (x,y,z).
// This will NOT produce valid results for any mesh other than the given
// sphere.
// Replace this with your own function for evaluating the implicit function
// values at the grid points using MLS
double wendlandWeight(double r, double h) {
    return pow(1 - r / h, 4) * (4 * r / h + 1);
}

void evaluateImplicitFunc()
{
    // dimension of the grid
    auto bb_min = grid_points.colwise().minCoeff().eval();
    auto bb_max = grid_points.colwise().maxCoeff().eval();
    
    // wendland radius
    wendlandRadius = 0.25 * (bb_max - bb_min).norm();
    // wendlandRadius = 100;

    // scalar values of the grid points (the implicit function values)
    grid_values.resize(resolution * resolution * resolution);

    // evaluate sphere's signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);
                // coordinates of the point at (x, y, z)
                Eigen::RowVector3d pi = grid_points.row(index);
                double xi = pi(0);
                double yi = pi(1);
                double zi = pi(2);
                // spatial index of the point at (x, y, z)
                int xIdx = min((int)floor((x - minX) / dimUnitX), resolutionSpatialIndex - 2);
                int yIdx = min((int)floor((y - minY) / dimUnitY), resolutionSpatialIndex - 2);
                int zIdx = min((int)floor((z - minZ) / dimUnitZ), resolutionSpatialIndex - 2);

                // find all points within wendlandRadius
                std::vector<int> neighbours;
                for (int localX = max(0, xIdx - 1); localX <= min(resolutionSpatialIndex - 2, xIdx + 1); localX ++) {
                    for (int localY = max(0, yIdx - 1); localY <= min(resolutionSpatialIndex - 2, yIdx + 1); localY ++) {
                        for (int localZ = max(0, zIdx - 1); localZ <= min(resolutionSpatialIndex - 2, zIdx + 1); localZ ++) {
                            for (int j : gridToVertices[localX][localY][localZ]) {
                                // check if the point is within the radius
                                Eigen::RowVector3d pj = constrained_points.row(j);
                                double dist = (pi - pj).norm();
                                double wj = wendlandWeight(dist, wendlandRadius);
                                if (dist < wendlandRadius && wj > 0.001) {
                                    neighbours.push_back(j);
                                }
                            }
                        }
                    }
                }

                // for (int j : gridToVertices[xIdx][yIdx][zIdx]) {
                //     // check if the point is within the radius
                //     Eigen::RowVector3d pj = constrained_points.row(j);
                //     double dist = (pi - pj).norm();
                //     double wj = wendlandWeight(dist, wendlandRadius);
                //     if (wj > 0) {
                //         neighbours.push_back(j);
                //     }
                // }

                // if no neighbours, set value to a large one
                if (neighbours.size() == 0) {
                    grid_values[index] = 100;
                    continue;
                }

                // determine the number of coefficients
                int numParams = 1;
                if (polyDegree == 0) {
                    numParams = 1;
                } else if (polyDegree == 1) {
                    numParams = 4;
                } else if (polyDegree == 2) {
                    numParams = 10;
                }

                Eigen::MatrixXd W = Eigen::MatrixXd::Zero(neighbours.size(), neighbours.size());
                Eigen::MatrixXd B = Eigen::MatrixXd::Zero(neighbours.size(), numParams);
                Eigen::VectorXd f = Eigen::VectorXd::Zero(neighbours.size(), 1);
                // for all neighbouring points, add one constraint
                for (int k = 0; k < neighbours.size(); k ++) {
                    int j = neighbours[k];

                    // get neighbour coordinates
                    Eigen::RowVector3d pj = constrained_points.row(j);
                    double xj = pj[0];
                    double yj = pj[1];
                    double zj = pj[2];
                    double dist = (pi - pj).norm();
                    
                    // compute wendland weight wj
                    double wj = wendlandWeight(dist, wendlandRadius);
                    W(k, k) = wj;
                    // compute b(pj)
                    if (polyDegree == 0) {
                        B.row(k)[0] = 1;
                    } else if (polyDegree == 1) {
                        B.row(k)[0] = 1;
                        B.row(k)[1] =  xj;
                        B.row(k)[2] = yj;
                        B.row(k)[3] = zj;
                    } else if (polyDegree == 2) {
                        B.row(k)[0] = 1;
                        B.row(k)[1] =  xj;
                        B.row(k)[2] = yj;
                        B.row(k)[3] = zj;
                        B.row(k)[4] =  xj * yj;
                        B.row(k)[5] = xj * zj;
                        B.row(k)[6] = yj* zj;
                        B.row(k)[7] = xj * xj;
                        B.row(k)[8] = yj * yj;
                        B.row(k)[9] = zj * zj;
                    }
                    // compute implicit function value fj
                    f[k] = constrained_values[j];
                }

                // solve least-square problem
                Eigen::VectorXd c = (W * B).colPivHouseholderQr().solve(W * f);
                // evaluate implicit function value at pi
                Eigen::VectorXd bi = Eigen::VectorXd(10);
                if (polyDegree == 0) {
                    bi[0] = 1;
                } else if (polyDegree == 1) {
                    bi[0] = 1;
                    bi[1] =  xi;
                    bi[2] = yi;
                    bi[3] = zi;
                } else if (polyDegree == 2) {
                    bi[0] = 1;
                    bi[1] =  xi;
                    bi[2] = yi;
                    bi[3] = zi;
                    bi[4] =  xi * yi;
                    bi[5] = xi * zi;
                    bi[6] = yi * zi;
                    bi[7] = xi * xi;
                    bi[8] = yi * yi;
                    bi[9] = zi * zi;
                }
                grid_values[index] = (bi.transpose() * c)[0];
            }
        }
    }
}

void evaluateImplicitFunc_PolygonSoup()
{
    // Replace with your code here, for "key == '5'"
    evaluateImplicitFunc();
}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines()
{
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                int index = x + resolution * (y + resolution * z);
                if (x < resolution - 1)
                {
                    int index1 = (x + 1) + y * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolution - 1)
                {
                    int index1 = x + (y + 1) * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolution - 1)
                {
                    int index1 = x + y * resolution + (z + 1) * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

// Estimation of the normals via PCA.
void pcaNormal()
{
    NP = -N; // to be replaced with your code
}

void create_spatial_index() {
    // spatial index resolution
    resolutionSpatialIndex = (resolution - 1) / 4 + 1;

    // get boundaries of the point cloud
    auto bb_min = P.colwise().minCoeff().eval();
    auto bb_max = P.colwise().maxCoeff().eval();

    // initialize epsilon
    diagonal = (bb_max - bb_min).norm();
    global_epsilon = 0.01 * diagonal;

    // compute dimension of a single grid
    dimUnitX = (bb_max[0] - bb_min[0]) / (double) (resolutionSpatialIndex - 1);
    dimUnitY = (bb_max[1] - bb_min[1]) / (double) (resolutionSpatialIndex - 1);
    dimUnitZ = (bb_max[2] - bb_min[2]) / (double) (resolutionSpatialIndex - 1);

    // spatial index boundaries given resolution
    minX = bb_min[0];
    minY = bb_min[1];
    minZ = bb_min[2];
    maxX = bb_max[0];
    maxY = bb_max[1];
    maxZ = bb_max[2];

    // initialize grid_to_vertices map
    gridToVertices.resize((int)resolutionSpatialIndex - 1, std::vector<std::vector<std::vector<int>>>((int)resolutionSpatialIndex - 1, std::vector<std::vector<int>>((int)resolutionSpatialIndex - 1)));
    // index the vertices
    for (int i = 0; i < P.rows(); i ++) {
        double x = P(i, 0);
        double y = P(i, 1);
        double z = P(i, 2);
        int xIdx = min((int)floor((x - minX) / dimUnitX), resolutionSpatialIndex - 2);
        int yIdx = min((int)floor((y - minY) / dimUnitY), resolutionSpatialIndex - 2);
        int zIdx = min((int)floor((z - minZ) / dimUnitZ), resolutionSpatialIndex - 2);
        gridToVertices[xIdx][yIdx][zIdx].push_back(i);
    }
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers)
{
    if (key == '1')
    {
        // Show imported points
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        viewer.data().point_size = 7;
        viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 0));
    }

    if (key == '2')
    {
        // load_grid_boundaries();
        create_spatial_index();
        gettimeofday(&startTime, NULL);

        // Show all constraints
        viewer.data().clear();
        viewer.core().align_camera_center(P);

        // initialize constrained points
        constrained_points.resize(3 * P.rows(), 3);
        constrained_values.resize(3 * P.rows());
        Eigen::MatrixXd constrained_points_colors(3 * P.rows(), 3);
        constrained_points_colors.setZero();

        // Add your code for computing auxiliary constraint points here
        for (int i = 0; i < P.rows(); i ++) {
            // load data
            double x = P(i, 0);
            double y = P(i, 1);
            double z = P(i, 2);
            // add a corresponding constraint
            constrained_points.row(3*i) = P.row(i);
            constrained_values(3*i) = 0;

            // compute the spatial index of the vertex
            int xIdx = floor((x - minX) / dimUnitX);
            int yIdx = floor((y - minY) / dimUnitY);
            int zIdx = floor((z - minZ) / dimUnitZ);
            int xStartIdx = max(0, xIdx - 1);
            int xEndIdx = min(resolutionSpatialIndex - 2, xIdx + 1);
            int yStartIdx = max(0, yIdx - 1);
            int yEndIdx = min(resolutionSpatialIndex - 2, yIdx + 1);
            int zStartIdx = max(0, zIdx - 1);
            int zEndIdx = min(resolutionSpatialIndex - 2, zIdx + 1);

            double epsilon = global_epsilon;
            // compute +epsilon
            Eigen::RowVector3d epsilonPlus = P.row(i) + epsilon * N.row(i);

            // check that p is the closest point to +epsilon
            int closestPt = i;
            double minDist = (epsilon * N.row(i)).norm();
            while (true) {
                // with spatial index
                for (int l = xStartIdx; l <= xEndIdx; l ++) {
                    for (int m = yStartIdx; m <= yEndIdx; m ++) {
                        for (int n = zStartIdx; n <= zEndIdx; n ++) {
                            for (int j : gridToVertices[l][m][n]) {
                                double tempDist = (P.row(j) - epsilonPlus).norm();
                                if (tempDist < minDist) {
                                    closestPt = j;
                                }
                            }
                        }
                    }
                }
                // // without spatial index
                // for (j = 0; j < P.rows(); j ++) {
                //     double tempX = P(j, 0);
                //     double tempY = P(j, 1);
                //     double tempZ = P(j, 2);
                //     double tempDistSqr = pow(tempX - epsilonPlusX, 2) + pow(tempY - epsilonPlusY, 2) + pow(tempZ - epsilonPlusZ, 2);
                //     if (tempDistSqr < minDistSqr) {
                //         closestPt = j;
                //     }
                // }

                if (closestPt == i) {
                    break;
                } else {
                    // halve epsilon
                    epsilon = epsilon / 2.0;
                    // recompute +epsilon
                    epsilonPlus = P.row(i) + epsilon * N.row(i);
                    minDist = (epsilon * N.row(i)).norm();
                }
            }

            // add +epsilon point to constrained_points
            constrained_points.row(3*i+1) = epsilonPlus;
            constrained_values(3*i+1) = epsilon;
            constrained_points_colors(3*i+1, 0) = 255;
            
            epsilon = global_epsilon;
            // compute -epsilon
            Eigen::RowVector3d epsilonMinus = P.row(i) - epsilon * N.row(i);
            // check that p is the closest point to -epsilon
            closestPt = i;
            minDist = (epsilon * N.row(i)).norm();
            // cout << "before second while" << endl;
            while (true) {
                // with spatial index
                for (int l = xStartIdx; l <= xEndIdx; l ++) {
                    for (int m = yStartIdx; m <= yEndIdx; m ++) {
                        for (int n = zStartIdx; n <= zEndIdx; n ++) {
                            for (int j : gridToVertices[l][m][n]) {
                                double tempDist = (P.row(j) - epsilonMinus).norm();
                                if (tempDist < minDist) {
                                    closestPt = j;
                                }
                            }
                        }
                    }
                }
                // // without spatial index
                // for (j = 0; j < P.rows(); j ++) {
                //     double tempX = P(j, 0);
                //     double tempY = P(j, 1);
                //     double tempZ = P(j, 2);
                //     double tempDistSqr = pow(tempX - epsilonMinusX, 2) + pow(tempY - epsilonMinusY, 2) + pow(tempZ - epsilonMinusZ, 2);
                //     if (tempDistSqr < minDistSqr) {
                //         closestPt = j;
                //     }
                // }

                if (closestPt == i) {
                    break;
                } else {
                    // halve epsilon
                    epsilon = epsilon / 2.0;
                    // recompute -epsilon
                    epsilonMinus = P.row(i) - epsilon * N.row(i);
                    minDist = (epsilon * N.row(i)).norm();
                }
            }
            // add -epsilon to constrained_points
            constrained_points.row(3*i+2) = epsilonMinus;
            constrained_values(3*i+2) = -epsilon;
            constrained_points_colors(3*i+2, 1) = 255;
        }

        // Add code for displaying all points, as above
        viewer.data().point_size = 7;
        viewer.data().add_points(constrained_points, constrained_points_colors);

        gettimeofday(&endTime, NULL);
        cout << "running time: " << (long) endTime.tv_sec - startTime.tv_sec << " milliseconds" << endl;
    }

    if (key == '3')
    {
        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);

        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc();

        // Add code for displaying points and lines

        // You can use the following example:
        /*** begin: sphere example, replace (at least partially) with your code ***/

        // get grid lines
        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i)
        {
            double value = grid_values(i);
            if (value < 0)
            {
                grid_colors(i, 1) = 1;
            }
            else
            {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }

        cout << grid_values << endl;

        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
                                grid_lines.block(0, 3, grid_lines.rows(), 3),
                                Eigen::RowVector3d(0.8, 0.8, 0.8));
        /*** end: sphere example ***/
    }

    if (key == '4')
    {
        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0))
        {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolution, resolution, resolution, V, F);
        if (V.rows() == 0)
        {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }

        igl::per_face_normals(V, F, FN);
        viewer.data().set_mesh(V, F);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);
    }

    if (key == '5')
    {
        // Use the structure for key=='3' but replace the function evaluateImplicitFunc();
        // with a function performing the approximation of the implicit surface from polygon soup
        // Ref: Chen Shen, James F. O’Brien, and Jonathan Richard Shewchuk. Interpolating and approximating implicit surfaces from polygon soup.

        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);

        // Make grid
        createGrid();

        // Evaluate implicit function --> Function to be modified here
        evaluateImplicitFunc_PolygonSoup();

        // get grid lines
        getLines();

        // Display the reconstruction
        callback_key_down(viewer, '4', modifiers);
    }

    if (key == '6' || key == '7' || key == '8')
    {
        // Implement PCA Normal Estimation --> Function to be modified here
        pcaNormal();

        // To use the normals estimated via PCA instead of the input normals and then restaurate the input normals
        Eigen::MatrixXd N_tmp = N;
        N = NP;

        switch (key)
        {
        case '6':
            callback_key_down(viewer, '2', modifiers);
            break;
        case '7':
            callback_key_down(viewer, '3', modifiers);
            break;
        case '8':
            callback_key_down(viewer, '3', modifiers);
            callback_key_down(viewer, '4', modifiers);
            break;
        default:
            break;
        }

        // Restore input normals
        N = N_tmp;
    }

    return true;
}

bool callback_load_mesh(Viewer &viewer, string filename)
{
    igl::readOFF(filename, P, F, N);
    callback_key_down(viewer, '1', 0);
    return true;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage ex2_bin <mesh.off>" << endl;
        igl::readOFF("../data/sphere.off", P, F, N);
    }
    else
    {
        // Read points and normals
        igl::readOFF(argv[1], P, F, N);
    }

    Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    viewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Expose variable directly ...
            ImGui::InputInt("Resolution", &resolution, 0, 0);
            if (ImGui::Button("Reset Grid", ImVec2(-1, 0)))
            {
                std::cout << "ResetGrid\n";
                // Recreate the grid
                createGrid();
                // Switch view to show the grid
                callback_key_down(viewer, '3', 0);
            }

            // TODO: Add more parameters to tweak here...
        }
    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
