#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/local_basis.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>


/*** insert any necessary libigl headers here ***/
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/cat.h>
#include <igl/dijkstra.h>

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
MatrixXd V;

// face array, #F x3
MatrixXi F;

// UV coordinates, #V x2
MatrixXd UV;

bool showingUV = false;
bool freeBoundary = false;
double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;

void Redraw()
{
	viewer.data().clear();

	if (!showingUV)
	{
		viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(false);

    if(UV.size() != 0)
    {
      viewer.data().set_uv(TextureResolution*UV);
      viewer.data().show_texture = true;
    }
	}
	else
	{
		viewer.data().show_texture = false;
		viewer.data().set_mesh(UV, F);
	}
}

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
	return false;
}

static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}
static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}

void ConvertConstraintsToMatrixForm(VectorXi indices, MatrixXd positions, Eigen::SparseMatrix<double> &C, VectorXd &d)
{
	// Convert the list of fixed indices and their fixed positions to a linear system
	// Hint: The matrix C should contain only one non-zero element per row and d should contain the positions in the correct order.
	
	// scalar replacement
	int vsize = V.rows();
	int wsize = indices.rows();

	// initialize C and d
	C.resize(2 * wsize, 2 * vsize);
	d.resize(2 * wsize, 1);

	// fill C
	for (int i = 0; i < wsize; i++) {
		C.insert(i, indices[i]) = 1;
		C.insert(wsize + i, vsize + indices[i]) = 1;
	}

	// fill d
	for (int i = 0; i < wsize; i++) {
		d[i] = positions.row(i)[0];
		d[wsize + i] = positions.row(i)[1];
	}
}

void computeParameterization(int type) {
	VectorXi fixed_UV_indices;
	MatrixXd fixed_UV_positions;

	SparseMatrix<double> A;
	VectorXd b;
	SparseMatrix<double> C;
	VectorXd d;

	// find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices
	std::vector<int> Bd;
	igl::boundary_loop(F, Bd);		// return the vertices of the longest boundary loop

	if (!freeBoundary) {
		// the boundary vertices should be fixed to positions on the unit disc. Find these position and
		// save them in the #V x 2 matrix fixed_UV_position.
		fixed_UV_indices.resize(Bd.size(), 1);
		for (int i = 0; i < Bd.size(); i ++) {
			fixed_UV_indices[i] = Bd[i];
		}
		fixed_UV_positions.resize(Bd.size(), 2);
		igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
	} else {
		// fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.
		// so far: two extreme points of the boundary loop

		// find the longest boundary loop
		VectorXi boundary_indices;
		MatrixXd boundary_positions;

		boundary_indices.resize(Bd.size(), 1);
		for (int i = 0; i < Bd.size(); i ++) {
			boundary_indices[i] = Bd[i];
		}
		boundary_positions.resize(Bd.size(), 2);
		igl::map_vertices_to_circle(V, boundary_indices, boundary_positions);

		// initialize
		fixed_UV_indices.resize(2, 1);
		fixed_UV_positions.resize(2, 2);

		// find the two vertices with the "longest" topological distance
		// consider all points
		set<int> targets;
		vector<vector<int>> VV;
		VectorXd min_distance;
		vector<double> max_distance;
		vector<int> max_distance_idx;
		VectorXi previous;
		igl::adjacency_list(F, VV);

		for (int source = 0; source < V.rows(); source ++) {
			igl::dijkstra(source, targets, VV, min_distance, previous);
			int idx;
			max_distance.push_back(min_distance.maxCoeff(&idx));
			max_distance_idx.push_back(idx);
		}

		double maxDist = 0;
		int maxIdx;
		int s, t;
		for (int i = 0; i < V.rows(); i ++) {
			if (max_distance[i] > maxDist) {
				maxDist = max_distance[i];
				maxIdx = max_distance_idx[i];
				s = i;
				t = maxIdx;
			}
		}

		fixed_UV_indices << s, t;
		fixed_UV_positions << -1, 0, 1, 0;
	}

	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);

	// Find the linear system for the parameterization (1- Tutte, 2- Harmonic, 3- LSCM, 4- ARAP)
	// and put it in the matrix A.
	// The dimensions of A should be 2#V x 2#V.
	if (type == '1') {
		// Add your code for computing uniform Laplacian for Tutte parameterization
		// Hint: use the adjacency matrix of the mesh

		A.resize(2*V.rows(), 2*V.rows());
		b.setZero(2 * V.rows());

		SparseMatrix<int> M;		// adjacency matrix of the mesh
		igl::adjacency_matrix(F, M);
		SparseVector<int> S;		// sum of each row of the adjacency matrix
		igl::sum(M, 1, S);		// dim = 1: sum along each row
		SparseMatrix<int> D;		// matrix with the elements of S on the diagonal
		igl::diag(S, D);

		SparseMatrix<int> L = M - D;

		// construct A by putting two L's on the diagonal
		for (int i = 0; i < L.rows(); i ++) {
			for (int j = 0; j < L.row(i).size(); j ++) {
				A.insert(i, j) = L.coeff(i, j);
				A.insert(L.rows() + i, L.rows() + j) = L.coeff(i, j);
			}
		}
	}

	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~

		A.resize(2*V.rows(), 2*V.rows());
		b.setZero(2 * V.rows());

		SparseMatrix<double> L;		// cotangent matrix of the mesh
		igl::cotmatrix(V, F, L);

		// construct A by putting two L's on the diagonal
		for (int i = 0; i < L.rows(); i ++) {
			for (int j = 0; j < L.row(i).size(); j ++) {
				A.insert(i, j) = L.coeff(i, j);
				A.insert(L.rows() + i, L.rows() + j) = L.coeff(i, j);
			}
		}
	}

	if (type == '3') {
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!

		A.resize(2*V.rows(), 2*V.rows());
		b.setZero(2 * V.rows());

		// compute surface gradients
		SparseMatrix<double> Dx, Dy;
		computeSurfaceGradientMatrix(Dx, Dy);
		SparseMatrix<double> DxT = Dx.transpose();
		SparseMatrix<double> DyT = Dy.transpose();

		// compute 2*area for the triangles
		VectorXd DAreaVec;
		igl::doublearea(V, F, DAreaVec);

		// compute area weight matrix
		Eigen::SparseMatrix<double> DAreaMat;
		DAreaMat.resize(F.rows(), F.rows());
		for (int i = 0; i < F.rows(); i ++) {
			DAreaMat.insert(i, i) = DAreaVec[i];
		}

		// construct linear system
		SparseMatrix<double> A11, A12, A21, A22;
		A11 = DxT*DAreaMat*Dx + DyT*DAreaMat*Dy;
		A12 = -DxT*DAreaMat*Dy + DyT*DAreaMat*Dx;
		A21 = -DyT*DAreaMat*Dx + DxT*DAreaMat*Dy;
		A22 = DxT*DAreaMat*Dx + DyT*DAreaMat*Dy;
		SparseMatrix<double> A_1, A_2;
		igl::cat(1, A11, A21, A_1);
		igl::cat(1, A12, A22, A_2);
		igl::cat(2, A_1, A_2, A);
	}

	if (type == '4') {
		// Add your code for computing ARAP system and right-hand side
		// Implement a function that computes the local step first
		// Then construct the matrix with the given rotation matrices

		A.resize(2*V.rows(), 2*V.rows());
		b.setZero(2 * V.rows());

		// initialize (press '3')
		// compute Jacobians (Dx * u, Dy * u; Dx * v, Dy * v)
		SparseMatrix<double> Dx, Dy;
		computeSurfaceGradientMatrix(Dx, Dy);
		SparseMatrix<double> DxT = Dx.transpose();
		SparseMatrix<double> DyT = Dy.transpose();
		VectorXd D11, D12, D21, D22;
		D11 = Dx*UV.col(0);		// Dx * u
		D12 = Dy*UV.col(0);		// Dy * u
		D21 = Dx*UV.col(1);		// Dx * v
		D22 = Dy*UV.col(1);		// Dy * v

		// compute closest rotation for each face
		MatrixXd R;
		R.resize(F.rows(), 4);
		for (int i = 0; i < F.rows(); i ++) {
			Eigen::Matrix2d J;		// Jacobian of face i
			J.resize(2, 2);
			J << D11[i], D12[i], D21[i], D22[i];

			Eigen::Matrix2d U, S, V;	// USV.T decomposition
			SSVD2x2(J, U, S, V);

			Eigen::Matrix2d tempR;		// closest rotation of face i
			tempR.resize(2, 2);
			Eigen::Matrix2d UVt; 
			UVt = U*V.transpose();
			if (UVt.determinant() >= 0) {
				tempR = UVt;
			} else {
				Eigen::Matrix2d temp;
				temp.resize(2, 2);
				temp << 1, 0, 0, -1;
				tempR = U*temp*V.transpose();
			}
			R.row(i) << tempR(0, 0), tempR(0, 1), tempR(1, 0), tempR(1, 1);
		}

		// construct linear system (copy from '3')
		VectorXd DAreaVec;
		igl::doublearea(V, F, DAreaVec);

		Eigen::SparseMatrix<double> AreaMat;
		AreaMat.resize(F.rows(), F.rows());
		for (int i = 0; i < F.rows(); i ++) {
			AreaMat.insert(i, i) = (DAreaVec[i] / 2.0);
		}

		SparseMatrix<double> A11, A12, A21, A22;
		A11 = DxT*AreaMat*Dx + DyT*AreaMat*Dy;
		A22 = DxT*AreaMat*Dx + DyT*AreaMat*Dy;
		A12.resize(V.rows(), V.rows());
		A21.resize(V.rows(), V.rows());
		SparseMatrix<double> A_1, A_2;
		igl::cat(1, A11, A21, A_1);
		igl::cat(1, A12, A22, A_2);
		igl::cat(2, A_1, A_2, A);

		VectorXd b1, b2;
		b1 = DxT*AreaMat*R.col(0) + DyT*AreaMat*R.col(1);
		b2 = DxT*AreaMat*R.col(2) + DyT*AreaMat*R.col(3);
		b.resize(2*V.rows());
		igl::cat(1, b1, b2, b);
	}

	// Solve the linear system.
	// Construct the system as discussed in class and the assignment sheet
	// Use igl::cat to concatenate matrices
	// Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail
	Eigen::SparseMatrix<double> ML, MR;
	Eigen::SparseMatrix<double, ColMajor> Zero, Mat;
	Zero.resize(C.rows(), C.rows());
	igl::cat(1, A, C, ML);
	Eigen::SparseMatrix<double> CT = C.transpose();
	igl::cat(1, CT, Zero, MR);
	igl::cat(2, ML, MR, Mat);

	SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solver;
	VectorXd bd;
	bd.resize(b.rows() + d.rows(), 1);
	bd << b, d;
	solver.analyzePattern(Mat);
	solver.factorize(Mat);
	VectorXd uvlamda;
	uvlamda = solver.solve(bd);

	// The solver will output a vector
	UV.resize(V.rows(), 2);
	UV.col(0) = uvlamda.block(0, 0, V.rows(), 1);
	UV.col(1) = uvlamda.block(V.rows(), 0, V.rows(), 1);
}

void colorCodeDistortion() {
	
}

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
	switch (key) {
	case '1':
	case '2':
	case '3':
	case '4':
		computeParameterization(key);
		break;
	case '5':
		// Add your code for detecting and displaying flipped triangles in the
		// UV domain here
		visualizeDistortion();
		break;
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
    if(showingUV)
    {
      temp2D = viewer.core();
      viewer.core() = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        temp3D = viewer.core();
        viewer.core() = temp2D;
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	Redraw();
	return true;
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  Redraw();
  viewer.core().align_camera_center(V);
  showingUV = false;

  return true;
}

bool callback_init(Viewer &viewer)
{
	temp3D = viewer.core();
	temp2D = viewer.core();
	temp2D.orthographic = true;

	return false;
}

int main(int argc,char *argv[]) {
  if(argc != 2) {
    cout << "Usage ex4_bin <mesh.off/obj>" << endl;
    load_mesh("../data/cathead.obj");
  }
  else
  {
    // Read points and normals
    load_mesh(argv[1]);
  }

	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Parmaterization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::Checkbox("Free boundary", &freeBoundary);

			// TODO: Add more parameters to tweak here...
		}
	};

  viewer.callback_key_pressed = callback_key_pressed;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_init = callback_init;

  viewer.launch();
}
