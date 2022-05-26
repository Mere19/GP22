#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>
#include <igl/invert_diag.h>
#include <igl/grad.h>
#include <igl/doublearea.h>
#include <igl/diag.h>
#include <igl/per_vertex_normals.h>

#include "Lasso.h"
#include "Colors.h"

//activate this for alternate UI (easier to debug)
//#define UPDATE_ONLY_ON_UP

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

//vertex array, #V x3
Eigen::MatrixXd V(0,3), V_cp(0, 3);
//face array, #F x3
Eigen::MatrixXi F(0,3);

//mouse interaction
enum MouseMode { SELECT, TRANSLATE, ROTATE, NONE };
MouseMode mouse_mode = NONE;
bool doit = false;
int down_mouse_x = -1, down_mouse_y = -1;

//for selecting vertices
std::unique_ptr<Lasso> lasso;
//list of currently selected vertices
Eigen::VectorXi selected_v(0,1);

//for saving constrained vertices
//vertex-to-handle index, #V x1 (-1 if vertex is free)
Eigen::VectorXi handle_id(0,1);
//list of all vertices belonging to handles, #HV x1
Eigen::VectorXi handle_vertices(0,1);
//centroids of handle regions, #H x1
Eigen::MatrixXd handle_centroids(0,3);
//updated positions of handle vertices, #HV x3
Eigen::MatrixXd handle_vertex_positions(0,3);
//index of handle being moved
int moving_handle = -1;
//rotation and translation for the handle being moved
Eigen::Vector3f translation(0,0,0);
Eigen::Vector4f rotation(0,0,0,1.);
typedef Eigen::Triplet<double> T;
//per vertex color array, #V x3
Eigen::MatrixXd vertex_colors;

//When false, use standard displacement vectors for details, when true use Deformation Transfer from part 2
bool use_deformation_transfer = false;

//function declarations (see below for implementation)
bool solve(Viewer& viewer);
void get_new_handle_locations();
Eigen::Vector3f computeTranslation (Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
Eigen::Vector4f computeRotation(Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
void compute_handle_centroids();
Eigen::MatrixXd readMatrix(const char *filename);

bool callback_mouse_down(Viewer& viewer, int button, int modifier);
bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y);
bool callback_mouse_up(Viewer& viewer, int button, int modifier);
bool callback_pre_draw(Viewer& viewer);
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);
void onNewHandleID();
void applySelection();

// for visualization
// int mesh = 1;
// enum MeshMode { SMOOTH, DEFORMED_SMOOTH, DEFORMED, ORIGINAL };
// MeshMode mesh_mode = ORIGINAL;
// int mesh = static_cast<int>(mesh_mode);
int mesh = 2;
bool default_algo = true;
// const char * items[]{"smooth", "deformed smooth", "deformed", };

Eigen::VectorXi free_vertices;
Eigen::MatrixXd Vs, Vds, Vd;
Eigen::SparseMatrix<double> Lw, M, Minv, Aff, Afc, nAff, nAfc;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> solver;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> non_default_solver;

void compute_free_vertices() {
  int c = 0;
  free_vertices.resize(V.rows() - handle_vertices.size());
  for (int i = 0; i < V.rows(); i ++) {
    if (handle_id[i] == -1) {
      free_vertices[c ++] = i;
    }
  }
}

void prefactor() {
  igl::cotmatrix(V, F, Lw);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
  igl::invert_diag(M, Minv);
  Eigen::SparseMatrix<double> A = Lw * Minv * Lw;    // compute A = Lw * Minv * Lw

  igl::slice(A, free_vertices, free_vertices, Aff);
  igl::slice(A, free_vertices, handle_vertices, Afc);

  solver.compute(Aff);    // left hand side
}

void remove_high_frequency_details() {
  // 2: remove high-frequency details
  Vs = V;   // smoothed mesh
  Eigen::MatrixXd Vc = igl::slice(V, handle_vertices, 1);
  Eigen::MatrixXd Vf_smoothed = solver.solve(-Afc * Vc);    // right hand side
  igl::slice_into(Vf_smoothed, free_vertices, 1, Vs);
}

Eigen::MatrixXd Nv, Ne, Nc, D;
Eigen::VectorXi n_list;
void compute_displacements() {
  // 4: compute displacement from V to Vs
  // compute local basis for B
  // Eigen::MatrixXd Nv;    // (unit?) vertex normal
  igl::per_vertex_normals(Vs, F, Nv);
  // Eigen::MatrixXd Ne;   // longest projection of incident edges on the tangent plane
  Ne.resize(Vs.rows(), 3);
  // Eigen::MatrixXd Nc;   // cross product of Nv and Nc
  Nc.resize(Vs.rows(), 3);
  
  // Eigen::VectorXi n_list;   // store the neighbour of the longest edge
  n_list.resize(V.rows());

  std::vector<std::vector<int>> adj_list;   // find the longest edge
  igl::adjacency_list(F, adj_list);
  for (int i = 0; i < V.rows(); i ++) {
    double maxNorm = -1;
    Eigen::RowVector3d nv = Nv.row(i);
    Eigen::RowVector3d ne;
    for (int j : adj_list[i]) {
      // Eigen::RowVector3d outgoingEdge = Vs.row(j) - Vs.row(i);
      // Eigen::RowVector3d projectionEdge = nv.dot(outgoingEdge)/(nv.squaredNorm())*nv;
      // Eigen::RowVector3d tempNe = outgoingEdge - projectionEdge;
      Eigen::RowVector3d tempNe = (Vs.row(j) - (nv.dot(Vs.row(j) - Vs.row(i)) * nv)) - Vs.row(i);
      if (tempNe.norm() > maxNorm) {
        maxNorm = tempNe.norm();
        tempNe.normalize();
        ne = tempNe;
        n_list[i] = j;
      }
    }
    Ne.row(i) = ne;
    Nc.row(i) = nv.cross(ne);
  }

  // Eigen::MatrixXd D;    // displacement vector in local basis
  D.resize(V.rows(), 3);
  for (int i = 0; i < V.rows(); i ++) {   // compute displacement and change basis
    Eigen::RowVector3d tempDisplacement = V.row(i) - Vs.row(i);
    // Eigen::Matrix3d fromLocal;
    Eigen::Matrix3d toLocal;
    // fromLocal.row(0) = Nv.row(i);
    // fromLocal.row(1) = Ne.row(i);
    // fromLocal.row(2) = Nc.row(i);
    // toLocal.row(0) = Nv.row(i);
    // toLocal.row(1) = Ne.row(i);
    // toLocal.row(2) = Nc.row(i);
    // fromLocal << Nv.row(i), Ne.row(i), Nc.row(i);
    // Eigen::Matrix3d fromLocalTranspose = fromLocal.transpose();    // [Nv.T Ne.T Nc.T]
    // cout << fromLocal << endl;
    // Eigen::Matrix3d fromWorld = fromLocal.inverse();
    // D.row(i) = (toLocal * tempDisplacement.transpose()).transpose();
    D(i, 0) = tempDisplacement.dot(Nv.row(i));
    D(i, 1) = tempDisplacement.dot(Ne.row(i));
    D(i, 2) = tempDisplacement.dot(Nc.row(i));
  }
}

void transfer_high_frequency_details() {
  // 4: transfer high frequency details
  // Vd = Vds;
  Vd.resize(V.rows(), 3);
  Nv.setZero();
  igl::per_vertex_normals(Vds, F, Nv);   // compute local basis for B'
  for (int i = 0; i < V.rows(); i ++) {
    Eigen::RowVector3d nv = Nv.row(i);
    int j = n_list[i];
    // Eigen::RowVector3d outgoingEdge = Vds.row(j) - Vds.row(i);
    // Eigen::RowVector3d projectionEdge = nv.dot(outgoingEdge)/(nv.squaredNorm())*nv;
    // Eigen::RowVector3d ne = outgoingEdge - projectionEdge;
    Eigen::RowVector3d ne = (Vds.row(j) - (nv.dot(Vds.row(j) - Vds.row(i)) * nv)) - Vds.row(i);
    ne.normalize();
    Eigen::RowVector3d nc = nv.cross(ne);

    // Eigen::Matrix3d fromLocal;    // add displacements back to B'
    // Eigen::Matrix3d toLocal;
    // fromLocal.row(0) = nv;
    // fromLocal.row(1) = ne;
    // fromLocal.row(2) = nc;
    // toLocal.row(0) = nv;
    // toLocal.row(1) = ne;
    // toLocal.row(2) = nc;
    // Eigen::Matrix3d toWorld = toLocal.inverse();
    // Eigen::Matrix3d fromLocalTranspose = fromLocal.transpose();    // [Nv.T Ne.T Nc.T]
    // Eigen::RowVector3d d = (toWorld * D.row(i).transpose()).transpose();
    Eigen::RowVector3d d = D(i, 0) * nv + D(i, 1) * ne + D(i, 2) * nc;
    Vd.row(i) = Vds.row(i) + d;
    // cout << (Vd.row(i) == V.row(i)) << endl;
  }
}

Eigen::SparseMatrix<double> G;
void compute_gradient_operator() {
  // compute
  igl::grad(V, F, G);

  // permute
  int fnum = F.rows();
  Eigen::SparseMatrix<double> G_permuted;
  G_permuted.resize(3*fnum, V.rows());
  // for (int i = 0, j = 0; i < 3*fnum, j < fnum; i +=3, j ++) {
  //   G_permuted.row(i) = G.row(j);
  //   G_permuted.row(i + 1) = G.row(j + fnum);
  //   G_permuted.row(i + 2) = G.row(j + 2 * fnum);
  // }
  Eigen::VectorXd rows_permuted, cols_permuted;
  rows_permuted.resize(3*fnum);
  cols_permuted.resize(V.rows());
  // for (int i = 0, j = 0; i < 3*fnum, j < fnum; i +=3, j ++) {
  //   rows_permuted[i] = j;
  //   rows_permuted[i + 1] = j + fnum;
  //   rows_permuted[i + 2] = j + 2*fnum;
  // }
  for(int i = 0; i < F.rows(); i++){
    // area_vec_resize[3 * i] = area_vec_resize[3 * i + 1] = area_vec_resize[3 * i + 2] = area_vec[i];
    rows_permuted[3*i] = i;
    rows_permuted[3*i + 1] = i + F.rows();
    rows_permuted[3*i + 2] = i + 2 * F.rows();
  }
  for (int i = 0; i < V.rows(); i ++) {
    cols_permuted[i] = i;
  }


  igl::slice(G, rows_permuted, cols_permuted, G);

  return ;
}

Eigen::VectorXd dblA;
Eigen::SparseMatrix<double> Diag;
void compute_weighting_matrix() {
  igl::doublearea(V, F, dblA);
  Eigen::VectorXd sglA = dblA / 2.0;
  Eigen::VectorXd sglA_vec;
  sglA_vec.resize(3*F.rows());

  for (int i = 0, j = 0; i < 3*F.rows(), j < F.rows(); i += 3, j ++) {
    sglA_vec[i] = sglA[j];
    sglA_vec[i+1] = sglA[j];
    sglA_vec[i+2] = sglA[j];
  }

  igl::diag(sglA_vec, Diag);

  return;
}

Eigen::MatrixXd S;
Eigen::MatrixXd Ns, Nds;
void compute_source_deformation_gradient() {
  S.resize(3 * F.rows(), 3);
  igl::per_face_normals(Vs, F, Ns);
  igl::per_face_normals(Vds, F, Nds);
  for (int i = 0, j = 0; i < F.rows(), j < 3*F.rows(); i ++, j += 3) {   // for each face, compute the source deformation gradient S
    Eigen::Vector3d q1 = Vs.row(F(i, 0));
    Eigen::Vector3d q2 = Vs.row(F(i, 1));
    Eigen::Vector3d q3 = Vs.row(F(i, 2));
    Eigen::Vector3d q1_ = Vds.row(F(i, 0));
    Eigen::Vector3d q2_ = Vds.row(F(i, 1));
    Eigen::Vector3d q3_ = Vds.row(F(i, 2));

    Eigen::Matrix3d Qds_T;
    Qds_T.col(0) = q1_ - q3_;
    Qds_T.col(1) = q2_ - q3_;
    Qds_T.col(2) = Nds.row(i);

    Eigen::Matrix3d Qs_T;
    Qs_T.col(0) = q1 - q3;
    Qs_T.col(1) = q2 - q3;
    Qs_T.col(2) = Ns.row(i);

    Eigen::Matrix3d S_i = Qds_T * Qs_T.inverse();
    Eigen::Matrix3d S_i_T = S_i.transpose();
    S.row(j) = S_i_T.row(0);
    S.row(j + 1) =  S_i_T.row(1);
    S.row(j + 2) =  S_i_T.row(2);
  }
  // cout << "[Finished] compute source deformation " << endl;
}

void non_default_prefactor() {
  Eigen::SparseMatrix<double> A = G.transpose() * Diag * G;   // compute G.T * D * G

  igl::slice(A, free_vertices, free_vertices, nAff);
  igl::slice(A, free_vertices, handle_vertices, nAfc);

  non_default_solver.compute(nAff);    // left hand side
}

void transfer_deformation() {
  Eigen::MatrixXd RHS1 = G.transpose() * Diag * S;
  Eigen::MatrixXd RHS1_sliced = igl::slice(RHS1, free_vertices, 1);
  // cout << "RHS1 dimension: " << RHS1.rows() << " x " << RHS1.cols() << endl;
  // cout << "RHS1_sliced dimension: " << RHS1_sliced.rows() << " x " << RHS1_sliced.cols() << endl;
  Eigen::MatrixXd MMM = nAfc * handle_vertex_positions;
  // cout << "MMM dimension: " << MMM.rows() << " x " << MMM.cols() << endl;
  Eigen::MatrixXd RHS = RHS1_sliced - MMM;
  Eigen::MatrixXd Vf_deformed = non_default_solver.solve(RHS);
  // cout << "Vf_deformed dimension: " << Vf_deformed.rows() << " x " << Vf_deformed.cols() << endl;
  // cout << "free_vertices dimension: " << free_vertices.rows() << " x " << free_vertices.cols() << endl;
  Vd = Vds;
  igl::slice_into(Vf_deformed, free_vertices, 1, Vd);
  igl::slice_into(handle_vertex_positions, handle_vertices, 1, Vd);
}

bool solve(Viewer& viewer)
{
  /**** Add your code for computing the deformation from handle_vertex_positions and handle_vertices here (replace following line) ****/
  // igl::slice_into(handle_vertex_positions, handle_vertices, 1, V);

  // 3: deform the smooth mesh
  // cout << "deform the smooth mesh ..." << endl;
  Vds = Vs;    // deformed mesh
  Eigen::MatrixXd Vf_deformed = solver.solve(-Afc * handle_vertex_positions);
  igl::slice_into(handle_vertex_positions, handle_vertices, 1, Vds);
  igl::slice_into(Vf_deformed, free_vertices, 1, Vds);

  if (default_algo) {
    transfer_high_frequency_details();
  } else {
    compute_source_deformation_gradient();
    transfer_deformation();
  }

  if (mesh == 0) {    // smooth mesh
    V = Vs;
  } else if (mesh == 1) {   // deformed smooth mesh
    V = Vds;
  } else if (mesh == 2) {   // deformed mesh
    V = Vd;
  }

  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);
  viewer.data().set_normals(N);

  return true;
};

void get_new_handle_locations()
{
  int count = 0;
  for (long vi = 0; vi < V.rows(); ++vi)
    if (handle_id[vi] >= 0)
    {
      Eigen::RowVector3f goalPosition = V.row(vi).cast<float>();
      if (handle_id[vi] == moving_handle) {
        if (mouse_mode == TRANSLATE)
          goalPosition += translation;
        else if (mouse_mode == ROTATE) {
            Eigen::RowVector3f  goalPositionCopy = goalPosition;
            goalPosition -= handle_centroids.row(moving_handle).cast<float>();
            igl::rotate_by_quat(goalPosition.data(), rotation.data(), goalPositionCopy.data());
            goalPosition = goalPositionCopy;
            goalPosition += handle_centroids.row(moving_handle).cast<float>();
        }
      }
      handle_vertex_positions.row(count++) = goalPosition.cast<double>();
    }
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  viewer.data().clear();
  viewer.data().set_mesh(V, F);

  viewer.core().align_camera_center(V);
  V_cp = V;
  handle_id.setConstant(V.rows(), 1, -1);
  // Initialize selector
  lasso = std::unique_ptr<Lasso>(new Lasso(V, F, viewer));

  selected_v.resize(0,1);

  return true;
}

int main(int argc, char *argv[])
{
  if(argc != 2) {
    cout << "Usage assignment5 mesh.off>" << endl;
    load_mesh("../data/woody-lo.off");
  }
  else
  {
    load_mesh(argv[1]);
  }

  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  menu.callback_draw_viewer_menu = [&]()
  {
    // Draw parent menu content
    menu.draw_viewer_menu();

    // Add new group
    if (ImGui::CollapsingHeader("Deformation Controls", ImGuiTreeNodeFlags_DefaultOpen))
    {
          int mouse_mode_type = static_cast<int>(mouse_mode);

          if (ImGui::Combo("Mouse Mode", &mouse_mode_type, "SELECT\0TRANSLATE\0ROTATE\0NONE\0"))
          {
            mouse_mode = static_cast<MouseMode>(mouse_mode_type);
          }

          if (ImGui::Button("Clear Selection", ImVec2(-1,0)))
          {
            selected_v.resize(0,1);
          }

          if (ImGui::Button("Apply Selection", ImVec2(-1,0)))
          {
            applySelection();
          }

          if (ImGui::Button("Clear Constraints", ImVec2(-1,0)))
          {
            handle_id.setConstant(V.rows(),1,-1);
          }

          if (ImGui::Button("Original Mesh", ImVec2(-1,0)))
          {
            V = V_cp;
            Eigen::MatrixXd N;
            igl::per_vertex_normals(V, F, N);
            viewer.data().set_normals(N);
          }

          if (ImGui::Button("Smooth Mesh", ImVec2(-1,0)))
          {
            V = Vs;
            Eigen::MatrixXd N;
            igl::per_vertex_normals(V, F, N);
            viewer.data().set_normals(N);
          }

          if (ImGui::Button("Deformed Smooth Mesh", ImVec2(-1,0)))
          {
            V = Vds;
            Eigen::MatrixXd N;
            igl::per_vertex_normals(V, F, N);
            viewer.data().set_normals(N);
          }

          if (ImGui::Button("Deformed Mesh", ImVec2(-1,0)))
          {
            V = Vd;
            Eigen::MatrixXd N;
            igl::per_vertex_normals(V, F, N);
            viewer.data().set_normals(N);
          }
          // if(ImGui::Checkbox("Deformation Transfer", &use_deformation_transfer)){}
          // ImGui::Combo("Mesh", &mesh, "SMOOTH\0DEFORMED_SMOOTH\0DEFORMED\0ORIGINAL\0");
    }
  };

  viewer.callback_key_down = callback_key_down;
  viewer.callback_mouse_down = callback_mouse_down;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_mouse_up = callback_mouse_up;
  viewer.callback_pre_draw = callback_pre_draw;

  viewer.data().point_size = 10;
  viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
  viewer.launch();
}


bool callback_mouse_down(Viewer& viewer, int button, int modifier)
{
  if (button == (int) Viewer::MouseButton::Right)
    return false;

  down_mouse_x = viewer.current_mouse_x;
  down_mouse_y = viewer.current_mouse_y;

  if (mouse_mode == SELECT)
  {
    if (lasso->strokeAdd(viewer.current_mouse_x, viewer.current_mouse_y) >=0)
      doit = true;
    else
      lasso->strokeReset();
  }
  else if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
  {
    int vi = lasso->pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
    if(vi>=0 && handle_id[vi]>=0)  //if a region was found, mark it for translation/rotation
    {
      moving_handle = handle_id[vi];
      get_new_handle_locations();
      doit = true;
    }
  }
  return doit;
}

int counter = 0;
bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y)
{
  counter = (counter + 1) % 5;
  if (counter != 0) {
    return false;
  }
  if (!doit)
    return false;
  if (mouse_mode == SELECT)
  {
    lasso->strokeAdd(mouse_x, mouse_y);
    return true;
  }
  if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
  {
    if (mouse_mode == TRANSLATE) {
      translation = computeTranslation(viewer,
                                       mouse_x,
                                       down_mouse_x,
                                       mouse_y,
                                       down_mouse_y,
                                       handle_centroids.row(moving_handle));
    }
    else {
      rotation = computeRotation(viewer,
                                 mouse_x,
                                 down_mouse_x,
                                 mouse_y,
                                 down_mouse_y,
                                 handle_centroids.row(moving_handle));
    }
    get_new_handle_locations();
#ifndef UPDATE_ONLY_ON_UP
    solve(viewer);
    down_mouse_x = mouse_x;
    down_mouse_y = mouse_y;
#endif
    return true;

  }
  return false;
}

bool callback_mouse_up(Viewer& viewer, int button, int modifier)
{
  if (!doit)
    return false;
  doit = false;
  if (mouse_mode == SELECT)
  {
    selected_v.resize(0,1);
    lasso->strokeFinish(selected_v);
    return true;
  }

  if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
  {
#ifdef UPDATE_ONLY_ON_UP
    if(moving_handle>=0)
      solve(viewer);
#endif
    translation.setZero();
    rotation.setZero(); rotation[3] = 1.;
    moving_handle = -1;

    compute_handle_centroids();

    return true;
  }

  return false;
};

bool callback_pre_draw(Viewer& viewer)
{
  // initialize vertex colors
  vertex_colors = Eigen::MatrixXd::Constant(V.rows(),3,.9);

  // first, color constraints
  int num = handle_id.maxCoeff();
  if (num == 0)
    num = 1;
  for (int i = 0; i<V.rows(); ++i)
    if (handle_id[i]!=-1)
    {
      int r = handle_id[i] % MAXNUMREGIONS;
      vertex_colors.row(i) << regionColors[r][0], regionColors[r][1], regionColors[r][2];
    }
  // then, color selection
  for (int i = 0; i<selected_v.size(); ++i)
    vertex_colors.row(selected_v[i]) << 131./255, 131./255, 131./255.;

  viewer.data().set_colors(vertex_colors);
  viewer.data().V_material_specular.fill(0);
  viewer.data().V_material_specular.col(3).fill(1);
  viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_DIFFUSE | igl::opengl::MeshGL::DIRTY_SPECULAR;


  //clear points and lines
  viewer.data().set_points(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXd::Zero(0,3));
  viewer.data().set_edges(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXi::Zero(0,3), Eigen::MatrixXd::Zero(0,3));

  //draw the stroke of the selection
  for (unsigned int i = 0; i<lasso->strokePoints.size(); ++i)
  {
    viewer.data().add_points(lasso->strokePoints[i],Eigen::RowVector3d(0.4,0.4,0.4));
    if(i>1)
      viewer.data().add_edges(lasso->strokePoints[i-1], lasso->strokePoints[i], Eigen::RowVector3d(0.7,0.7,0.7));
  }

  // update the vertex position all the time
  viewer.data().V.resize(V.rows(),3);
  viewer.data().V << V;

  viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_POSITION;

#ifdef UPDATE_ONLY_ON_UP
  //draw only the moving parts with a white line
  if (moving_handle>=0)
  {
    Eigen::MatrixXd edges(3*F.rows(),6);
    int num_edges = 0;
    for (int fi = 0; fi<F.rows(); ++fi)
    {
      int firstPickedVertex = -1;
      for(int vi = 0; vi<3 ; ++vi)
        if (handle_id[F(fi,vi)] == moving_handle)
        {
          firstPickedVertex = vi;
          break;
        }
      if(firstPickedVertex==-1)
        continue;


      Eigen::Matrix3d points;
      for(int vi = 0; vi<3; ++vi)
      {
        int vertex_id = F(fi,vi);
        if (handle_id[vertex_id] == moving_handle)
        {
          int index = -1;
          // if face is already constrained, find index in the constraints
          (handle_vertices.array()-vertex_id).cwiseAbs().minCoeff(&index);
          points.row(vi) = handle_vertex_positions.row(index);
        }
        else
          points.row(vi) =  V.row(vertex_id);

      }
      edges.row(num_edges++) << points.row(0), points.row(1);
      edges.row(num_edges++) << points.row(1), points.row(2);
      edges.row(num_edges++) << points.row(2), points.row(0);
    }
    edges.conservativeResize(num_edges, Eigen::NoChange);
    viewer.data().add_edges(edges.leftCols(3), edges.rightCols(3), Eigen::RowVector3d(0.9,0.9,0.9));

  }
#endif
  return false;

}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers)
{
  bool handled = false;
  if (key == 'S')
  {
    mouse_mode = SELECT;
    handled = true;
  }

  if ((key == 'T') && (modifiers == IGL_MOD_ALT))
  {
    mouse_mode = TRANSLATE;
    handled = true;
  }

  if ((key == 'R') && (modifiers == IGL_MOD_ALT))
  {
    mouse_mode = ROTATE;
    handled = true;
  }
  if (key == 'A')
  {
    applySelection();
    callback_key_down(viewer, '1', 0);
    handled = true;
  }
  // if (key == '2')
  // viewer.core().align_camera_position(B);
  //  }

  //viewer.ngui->refresh();
  return handled;
}

void onNewHandleID()
{
  //store handle vertices too
  int numFree = (handle_id.array() == -1).cast<int>().sum();
  int num_handle_vertices = V.rows() - numFree;
  handle_vertices.setZero(num_handle_vertices);
  handle_vertex_positions.setZero(num_handle_vertices,3);

  int count = 0;
  for (long vi = 0; vi<V.rows(); ++vi)
    if(handle_id[vi] >=0)
      handle_vertices[count++] = vi;

  compute_handle_centroids();
  compute_free_vertices();
  prefactor();
  remove_high_frequency_details();

  if (default_algo) {
    compute_displacements();
  } else {
    // cout << "1" << endl;
    compute_gradient_operator();
    // cout << "2" << endl;
    compute_weighting_matrix();
    // cout << "3" << endl;
    non_default_prefactor();
    // cout << "4" << endl;
  }
}

void applySelection()
{
  int index = handle_id.maxCoeff()+1;
  for (int i =0; i<selected_v.rows(); ++i)
  {
    const int selected_vertex = selected_v[i];
    if (handle_id[selected_vertex] == -1)
      handle_id[selected_vertex] = index;
  }
  selected_v.resize(0,1);

  onNewHandleID();
}

void compute_handle_centroids()
{
  //compute centroids of handles
  int num_handles = handle_id.maxCoeff()+1;
  handle_centroids.setZero(num_handles,3);

  Eigen::VectorXi num; num.setZero(num_handles,1);
  for (long vi = 0; vi<V.rows(); ++vi)
  {
    int r = handle_id[vi];
    if ( r!= -1)
    {
      handle_centroids.row(r) += V.row(vi);
      num[r]++;
    }
  }

  for (long i = 0; i<num_handles; ++i)
    handle_centroids.row(i) = handle_centroids.row(i).array()/num[i];

}

//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector3f computeTranslation (Viewer& viewer,
                                    int mouse_x,
                                    int from_x,
                                    int mouse_y,
                                    int from_y,
                                    Eigen::RowVector3d pt3D)
{
  Eigen::Matrix4f modelview = viewer.core().view;// * viewer.data().model;
  //project the given point (typically the handle centroid) to get a screen space depth
  Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
                                      modelview,
                                      viewer.core().proj,
                                      viewer.core().viewport);
  float depth = proj[2];

  double x, y;
  Eigen::Vector3f pos1, pos0;

  //unproject from- and to- points
  x = mouse_x;
  y = viewer.core().viewport(3) - mouse_y;
  pos1 = igl::unproject(Eigen::Vector3f(x,y,depth),
                        modelview,
                        viewer.core().proj,
                        viewer.core().viewport);


  x = from_x;
  y = viewer.core().viewport(3) - from_y;
  pos0 = igl::unproject(Eigen::Vector3f(x,y,depth),
                        modelview,
                        viewer.core().proj,
                        viewer.core().viewport);

  //translation is the vector connecting the two
  Eigen::Vector3f translation = pos1 - pos0;
  return translation;

}


//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector4f computeRotation(Viewer& viewer,
                                int mouse_x,
                                int from_x,
                                int mouse_y,
                                int from_y,
                                Eigen::RowVector3d pt3D)
{

  Eigen::Vector4f rotation;
  rotation.setZero();
  rotation[3] = 1.;

  Eigen::Matrix4f modelview = viewer.core().view;// * viewer.data().model;

  //initialize a trackball around the handle that is being rotated
  //the trackball has (approximately) width w and height h
  double w = viewer.core().viewport[2]/8;
  double h = viewer.core().viewport[3]/8;

  //the mouse motion has to be expressed with respect to its center of mass
  //(i.e. it should approximately fall inside the region of the trackball)

  //project the given point on the handle(centroid)
  Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
                                      modelview,
                                      viewer.core().proj,
                                      viewer.core().viewport);
  proj[1] = viewer.core().viewport[3] - proj[1];

  //express the mouse points w.r.t the centroid
  from_x -= proj[0]; mouse_x -= proj[0];
  from_y -= proj[1]; mouse_y -= proj[1];

  //shift so that the range is from 0-w and 0-h respectively (similarly to a standard viewport)
  from_x += w/2; mouse_x += w/2;
  from_y += h/2; mouse_y += h/2;

  //get rotation from trackball
  Eigen::Vector4f drot = viewer.core().trackball_angle.coeffs();
  Eigen::Vector4f drot_conj;
  igl::quat_conjugate(drot.data(), drot_conj.data());
  igl::trackball(w, h, float(1.), rotation.data(), from_x, from_y, mouse_x, mouse_y, rotation.data());

  //account for the modelview rotation: prerotate by modelview (place model back to the original
  //unrotated frame), postrotate by inverse modelview
  Eigen::Vector4f out;
  igl::quat_mult(rotation.data(), drot.data(), out.data());
  igl::quat_mult(drot_conj.data(), out.data(), rotation.data());
  return rotation;
}