#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any libigl headers here ***/
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/facet_components.h>
#include <igl/jet.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

void subdivide_sqrt3(const Eigen::MatrixXd &V,
					 const Eigen::MatrixXi &F,
					 Eigen::MatrixXd &Vout,
					 Eigen::MatrixXi &Fout){
    // copy old vertices to Vout
    for (int i = 0; i < V.rows(); i ++) {
        Vout(i) = V(i);
    }
    std::cout << Vout << std::endl;
    // copy old faces to Fout
    for (int i = 0; i < F.rows(); i ++) {
        Fout(i) = F(i);
    }
    std::cout << Fout << std::endl;

    // locate new vertices
    int j = V.rows();
    int k = F.rows();
    for (int i = 0; i < F.rows(); i ++) {
        double x = (V(F(i, 0), 0) + V(F(i, 1), 0) + V(F(i, 2), 0)) / 3;
        double y = (V(F(i, 0), 1) + V(F(i, 1), 1) + V(F(i, 2), 1)) / 3;
        double z = (V(F(i, 0), 2) + V(F(i, 1), 2) + V(F(i, 2), 2)) / 3;

        // add new vertex to Vout
        Vout(j, 0) = x;
        Vout(j, 1) = y;
        Vout(j, 2) = z;
        std::cout << Vout(j) << std::endl;
        j ++;

        std::cout << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << std::endl;

        // add new faces to Fout
        Fout(k, 0) = F(i, 0);
        Fout(k, 1) = F(i, 1);
        Fout(k++, 2) = j - 1;
        std::cout << Fout(k - 1, 0) << Fout(k - 1, 1) << Fout(k - 1, 2) << std::endl;

        Fout(k, 0) = F(i, 1);
        Fout(k, 1) = F(i, 2);
        Fout(k++, 2) = j - 1;
        std::cout << Fout(k - 1, 0) << Fout(k - 1, 1) << Fout(k - 1, 2) << std::endl;

        Fout(k, 0) = F(i, 2);
        Fout(k, 1) = F(i, 0);
        Fout(k++, 2) = j - 1;
        std::cout << Fout(k - 1, 0) << Fout(k - 1, 1) << Fout(k - 1, 2) << std::endl;

    }
}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to face relations here;
        // store in VF,VFi.
        igl::vertex_triangle_adjacency(V, F, VF, VFi);
        for (int i = 0; i < VF.size(); i++) {
            std::cout << "| vertex " << i << ": ";
            for (int fidx : VF[i]) {
                std::cout << fidx << " ";
            }
            std::cout << "|" << std::endl;
        }
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to vertex relations here:
        // store in VV.
        igl::adjacency_list(F, VV);
        for (int i = 0; i < VV.size(); i++) {
            std::cout << "| vertex " << i << ": ";
            for (int vidx : VV[i]) {
                std::cout << vidx << " ";
            }
            std::cout << "|"<< std::endl;
        }
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        FN.setZero(F.rows(),3);
        // Add your code for computing per-face normals here: store in FN.
        igl::per_face_normals(V, F, FN);
        // Set the viewer normals.
        viewer.data().set_normals(FN);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-vertex normals here: store in VN.
        igl::per_vertex_normals(V, F, VN);
        // Set the viewer normals.
        viewer.data().set_normals(VN);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-corner normals here: store in CN.
        igl::per_corner_normals(V, F, 75, CN);
        //Set the viewer normals
        viewer.data().set_normals(CN);
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        component_colors_per_face.setZero(F.rows(),3);

        // Add your code for computing per-face connected components here:
        // store the component labels in cid.
        igl::facet_components(F, cid);
        int max = 0;
        for (int i = 0; i < cid.size(); i ++) {
            if (cid[i] > max) max = cid[i];
        }
        std::cout << "number of components: " << max + 1 << std::endl;
        std::vector<int> count_component(max + 1, 0);
        for (int i = 0; i < cid.size(); i ++) {
            count_component[cid[i]] ++;
        }
        for (int i = 0; i < count_component.size(); i ++) {
            std::cout << "component " << i << ": " << count_component[i] << std::endl;
        }

        // Compute colors for the faces based on components, storing them in
        // component_colors_per_face.
        igl::jet(cid, true, component_colors_per_face);
        // Set the viewer colors
        viewer.data().set_colors(component_colors_per_face);
    }

    if (key == '7') {
        // std::cout << "in 7" << std::endl;
		Eigen::MatrixXd Vout(V.rows() + F.rows(), 3);
		Eigen::MatrixXi Fout(3 * F.rows(), 3);
        // Eigen::MatrixXi Fout;
        // Fill the subdivide_sqrt3() function with your code for sqrt(3) subdivision.
		subdivide_sqrt3(V,F,Vout,Fout);
        // Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
    }

    return true;
}

bool load_mesh(Viewer& viewer,string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
  igl::readOFF(filename,V,F);
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
        filename = std::string(argv[1]); // Mesh provided as command line argument
    }
    else {
        filename = std::string("../data/bunny.off"); // Default mesh
    }
	
    load_mesh(viewer,filename,V,F);

    callback_key_down(viewer, '2', 0);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
    viewer.launch();
}
