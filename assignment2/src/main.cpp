#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <windows.h>
#include <string>
#include <filesystem>
#include <igl/per_face_normals.h>
#include <igl/floor.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/slice.h>

using namespace std;

using Viewer = igl::opengl::glfw::Viewer;

Eigen::RowVector3d BLACK_COLOR = { 0 ,0, 0 };
Eigen::RowVector3d GREEN_COLOR = { 1 ,0, 0 };
Eigen::RowVector3d RED_COLOR = { 0 ,1, 0 };
Eigen::RowVector3d BLUE_COLOR = { 0 ,0, 1 };


// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
int polyDegree = 0;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;

// Parameter: grid resolution
int resolutionX = 20 , resolutionY = 20 , resolutionZ = 20;

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

// Bounding box diag.
double diag, spatialRate = 0.1;

// Spatial index length of each dimension
int spatialX, spatialY, spatialZ;

// Constrained points list indices
std::vector<std::vector<int> > spatialPoints;

std::string mesh_name;




// Functions
void createGrid();
void evaluateImplicitFunc();
void getLines();
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);

// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid() {
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines. resize(0, 6);
    grid_values.resize(0);
    V. resize(0, 3);
    F. resize(0, 3);
    FN.resize(0, 3);

    // Grid bounds: axis-aligned bounding box
    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();

    // Bounding box dimensions
    Eigen::RowVector3d dim = bb_max - bb_min;

    // Grid spacing
    const double dx = dim[0] / (double)(resolutionX - 1);
    const double dy = dim[1] / (double)(resolutionY - 1);
    const double dz = dim[2] / (double)(resolutionZ - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolutionX * resolutionY * resolutionZ, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolutionX; ++x) {
        for (unsigned int y = 0; y < resolutionY; ++y) {
            for (unsigned int z = 0; z < resolutionZ; ++z) {
                // Linear index of the point at (x,y,z)
                int index = x + resolutionY * (y + resolutionZ * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }
}

void find_nearby_points(Eigen::RowVector3d point, double h, std::vector<int>& result_vec, std::vector<double>& d_vec) {
    result_vec.clear();
    d_vec.clear();
    double unit_length = spatialRate * diag;
    Eigen::RowVector3d p_dim = (point - P.colwise().minCoeff())/unit_length;
    int num_cells = ceil(h / unit_length);
    int x1 = max(0, int(p_dim[0]) - num_cells); int x2 = min(spatialX, int(p_dim[0]) + num_cells + 1);
    int y1 = max(0, int(p_dim[1]) - num_cells); int y2 = min(spatialY, int(p_dim[1]) + num_cells + 1);
    int z1 = max(0, int(p_dim[2]) - num_cells); int z2 = min(spatialZ, int(p_dim[2]) + num_cells + 1);
    

    
    for (int i = x1; i < x2; i++) {
        for (int j = y1; j < y2; j++) {
            for (int k = z1; k < z2; k++) {
                for (int it = 0; it < spatialPoints[i + spatialX * (j + spatialY * k)].size(); it++) {
                    int it_idx = spatialPoints[i + spatialX * (j + spatialY * k)][it];
                    double distance = (constrained_points.row(it_idx) - point).norm();
                    if (distance <= h) {
                        result_vec.push_back(it_idx);
                        d_vec.push_back(distance);
                    }
                    distance = (constrained_points.row(P.rows()+it_idx) - point).norm();
                    if (distance <= h) {
                        result_vec.push_back(P.rows()+it_idx);
                        d_vec.push_back(distance);
                    }
                    distance = (constrained_points.row(P.rows()*2+it_idx) - point).norm();
                    if (distance <= h) {
                        result_vec.push_back(P.rows()*2+it_idx);
                        d_vec.push_back(distance);
                    }
                }
            }
        }
    }
}

// Function for evaluating the implicit function
// values at the grid points using MLS
void evaluateImplicitFunc() {
    double localWendlandRadius = wendlandRadius * diag;
    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolutionX * resolutionY * resolutionZ);

    // Evaluate sphere's signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolutionX; ++x) {
        for (unsigned int y = 0; y < resolutionY; ++y) {
            for (unsigned int z = 0; z < resolutionZ; ++z) {
                // Linear index of the point at (x,y,z)
                int index = x + resolutionX * (y + resolutionY * z);
                std::vector<int> result_vec;
                std::vector<double> d_vec;
                find_nearby_points(grid_points.row(index), localWendlandRadius, result_vec, d_vec);
                if (result_vec.size() == 0) {
                    grid_values[index] = 100;
                }
                else {
                    Eigen::VectorXd r_vec = Eigen::VectorXd::Map(d_vec.data(), d_vec.size());
                    
                    Eigen::MatrixXd nearby_points;
                    Eigen::VectorXd fi;
                    Eigen::VectorXi col_array(3);
                    col_array << 0, 1, 2;
                    igl::slice(constrained_points, Eigen::VectorXi::Map(result_vec.data(), result_vec.size()), col_array, nearby_points);
                    col_array.resize(1);
                    col_array << 0;
                    igl::slice(constrained_values, Eigen::VectorXi::Map(result_vec.data(), result_vec.size()), col_array, fi);

                    Eigen::MatrixXd M_squares = nearby_points.cwiseProduct(nearby_points);
                    Eigen::MatrixXd nearby_points_r(nearby_points.rows(), 3);
                    nearby_points_r << nearby_points.col(1), nearby_points.col(2), nearby_points.col(0);
                    Eigen::MatrixXd M_products = nearby_points.cwiseProduct(nearby_points_r);

                    Eigen::MatrixXd A;
                    Eigen::VectorXd bx;
                    double x = grid_points(index, 0);
                    double y = grid_points(index, 1);
                    double z = grid_points(index, 2);
                    if (polyDegree == 0) {
                        A = Eigen::MatrixXd::Ones(nearby_points.rows(), 1);
                        bx.resize(1);
                        bx << 1;

                    }
                    else if (polyDegree == 1) {
                        A.resize(nearby_points.rows(), 4);
                        A << Eigen::MatrixXd::Ones(nearby_points.rows(), 1), nearby_points;
                        bx.resize(4);
                        bx << 1, x, y, z;
                    }
                    else {
                        A.resize(nearby_points.rows(), 10);
                        A << Eigen::MatrixXd::Ones(nearby_points.rows(), 1), nearby_points, M_squares, M_products;
                        bx.resize(10);
                        bx << 1, x, y, z, pow(x, 2), pow(y, 2), pow(z, 2), (x * y), (y * z), (z * x);
                    }

                    Eigen::VectorXd W = (1 - (r_vec.array() / localWendlandRadius)).pow(4) * (4 * r_vec.array() / localWendlandRadius + 1);
                    Eigen::VectorXd c = (A.transpose() * W.asDiagonal() * A).ldlt().solve(A.transpose() * W.asDiagonal() * fi);
                    grid_values[index] = bx.dot(c);
                }
            }
        }
    }
}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines() {
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x<resolutionX; ++x) {
        for (unsigned int y = 0; y < resolutionY; ++y) {
            for (unsigned int z = 0; z < resolutionZ; ++z) {
                int index = x + resolutionX * (y + resolutionY * z);
                if (x < resolutionX - 1) {
                    int index1 = (x + 1) + y * resolutionX + z * resolutionY * resolutionZ;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolutionY - 1) {
                    int index1 = x + (y + 1) * resolutionX + z * resolutionY * resolutionZ;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolutionZ - 1) {
                    int index1 = x + y * resolutionX + (z + 1) * resolutionY * resolutionZ;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

int findPoint(Eigen::RowVectorXd curr) {
    int min_index = -1;
    double d;
    double d_min = diag;
    double unit_length = spatialRate * d_min;
    Eigen::RowVector3d p_dim = (curr - P.colwise().minCoeff()) / unit_length;
    for (int i = max(0, int(p_dim[0]) - 1); i < min(spatialX, int(p_dim[0]) + 2); i++) {
        for (int j = max(0, int(p_dim[1]) - 1); j < min(spatialY, int(p_dim[1]) + 2); j++) {
            for (int k = max(0, int(p_dim[2]) - 1); k < min(spatialZ, int(p_dim[2]) + 2); k++) {
                std::vector<int> spatialPointsVal = spatialPoints[i + spatialX * (j + spatialY * k)];
                for (int it = 0; it < spatialPointsVal.size(); it++) {
                    int it_idx = spatialPointsVal[it];
                    d = (P.row(it_idx) - curr).norm();
                    if (d < d_min) {
                        d_min = d;
                        min_index = it_idx;
                    }
                }
            }
        }
    }
    return min_index;
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        // Show imported points
        viewer.data().clear();
        viewer.core.align_camera_center(P);
        viewer.data().point_size = 5;
        viewer.data().add_points(P, BLACK_COLOR);
    }

    if (key == '2') {
        // Show all constraints
        viewer.data().clear();
        viewer.core.align_camera_center(P);

        diag = igl::bounding_box_diagonal(P);

        // Implementing a spatial index to accelerate neighbour calculations
        spatialPoints.clear();
        double unit_length = spatialRate * diag;
        Eigen::VectorwiseOp<Eigen::MatrixXd, 0> cols = P.colwise();
        Eigen::RowVector3d dim_length = (cols.maxCoeff() - cols.minCoeff()) / unit_length;

        spatialX = int(dim_length[0]) + 1; spatialY = int(dim_length[1]) + 1; spatialZ = int(dim_length[2]) + 1;
        Eigen::MatrixXd P_idx = (P - cols.minCoeff().replicate(P.rows(), 1)) / unit_length;

        spatialPoints.resize(spatialX * spatialY * spatialZ);
        for (int i = 0; i < P_idx.rows(); i++) {
            int spatial_idx = int(P_idx(i, 0)) + spatialX * int((P_idx(i, 1)) + spatialY * int(P_idx(i, 2)));
            spatialPoints[spatial_idx].push_back(i);
        }

        constrained_points.resize(P.rows() * 3, 3);
        constrained_values.setZero(P.rows() * 3);
        viewer.data().point_size = 5;
        viewer.data().add_points(P, BLUE_COLOR);
       
        double epsilon;
        Eigen::RowVectorXd positiveConstranined(3);
        Eigen::RowVectorXd negativeConstranined(3);
        for (int i = 0; i < P.rows(); i++) {
            epsilon = 0.03 * diag;
            positiveConstranined = P.row(i) + epsilon * N.row(i).normalized();
            while (findPoint(positiveConstranined) != i) {
                epsilon /= 2;
                positiveConstranined = P.row(i) + epsilon * N.row(i).normalized();
            }
            constrained_points.row(P.rows() + i) = positiveConstranined;
            constrained_values(P.rows() + i) = epsilon;
        }
        viewer.data().add_points(constrained_points.block(P.rows(), 0, P.rows(), 3), GREEN_COLOR);

        int P_multiply_idx = 2 * P.rows();
        for (int i = 0; i < P.rows(); i++) {
            epsilon = 0.01 * diag;
            negativeConstranined = P.row(i) - epsilon * N.row(i).normalized();
            while (findPoint(negativeConstranined) != i) {
                epsilon /= 2;
                negativeConstranined = P.row(i) - epsilon * N.row(i).normalized();
            }
            constrained_points.row(P_multiply_idx + i) = negativeConstranined;
            constrained_values(P_multiply_idx + i) = (-1) * epsilon;
        }
        viewer.data().add_points(constrained_points.block(P.rows() * 2, 0, P.rows(), 3), RED_COLOR);
        
    }

    if (key == '3') {
        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core.align_camera_center(P);

        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc();

        // get grid lines
        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i) {
            double value = grid_values(i);
            if (value < 0) {
                grid_colors(i, 1) = 1;
            }
            else {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }

        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
                              grid_lines.block(0, 3, grid_lines.rows(), 3),
                              Eigen::RowVector3d(0.8, 0.8, 0.8));
    }

    if (key == '4') {
        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0)) {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolutionX, resolutionY, resolutionZ, V, F);
        if (V.rows() == 0) {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }

        igl::per_face_normals(V, F, FN);
        viewer.data().set_mesh(V, F);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);
        
        if (resolutionX == resolutionY && resolutionX == resolutionZ) {
            igl::writeOFF("F:/Repos/geometryprocessing2021-sagigvili/assignment2/results/" + mesh_name + "_" + to_string(resolutionX) + "_" + to_string(wendlandRadius) + "_" + to_string(polyDegree) +".off", V, F);
        }
        else {
            igl::writeOFF("F:/Repos/geometryprocessing2021-sagigvili/assignment2/results/" + mesh_name + "_" + to_string(resolutionX) + "_" + to_string(resolutionY) + "_" + to_string(resolutionZ) + "_" + to_string(wendlandRadius) + "_" + to_string(polyDegree) + ".off", V, F);
        }
        
    }

    return true;
}

bool callback_load_mesh(Viewer& viewer,string filename)
{
  igl::readOFF(filename,P,F,N);
  std::size_t found = filename.find_last_of("/\\");
  mesh_name = filename.substr(found + 1);
  stringstream ss(mesh_name);
  std::string token;
  getline(ss, token, '.');
  mesh_name = token;
  cout << mesh_name << endl;
  callback_key_down(viewer,'1',0);
  return true;
}

std::string workingdir()
{
    char buf[256];
    GetCurrentDirectoryA(256, buf);
    return std::string(buf) + '\\';
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage ex2_bin <mesh.off>" << endl;
        igl::readOFF("../../../data/cat.off",P,F,N);
        mesh_name = "cat";
    }
	  else
	  {
		  // Read points and normals
		  igl::readOFF(argv[1],P,F,N);
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
        ImGui::InputInt("Resolution X", &resolutionX, 0, 0);
        ImGui::InputInt("Resolution Y", &resolutionY, 0, 0);
        ImGui::InputInt("Resolution Z", &resolutionZ, 0, 0);
        ImGui::InputDouble("wendlandRadiusRate", &wendlandRadius, 0, 0);
        ImGui::InputDouble("SpatialRate", &spatialRate, 0, 0);
        ImGui::InputInt("polyDegree", &polyDegree, 0, 0);
        if (ImGui::Button("Reset Grid", ImVec2(-1,0)))
        {
          std::cout << "ResetGrid\n";
          // Recreate the grid
          createGrid();
          // Switch view to show the grid
          callback_key_down(viewer,'3',0);
        }

        if (ImGui::Button("Upload Cloud", ImVec2(-1, 0)))
        {
            std::string fname = igl::file_dialog_open();
            std::cout << fname << " cloud has been uploaded\n";
            callback_load_mesh(viewer, fname);
        }

        // TODO: Add more parameters to tweak here...
      }

    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
