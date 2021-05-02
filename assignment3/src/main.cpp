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
#include <igl/adjacency_list.h>

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

// #F x 1 Vector - Distortion measurement for each face 
Eigen::VectorXd distortion;

// #F x 3 Matrix - distortion colors
Eigen::MatrixXd distortionMap;

bool showingUV = false;
bool freeBoundary = true;
bool conformal = false;
bool isometric = false;
int distortionType = -1;
double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;

void Redraw()
{
    viewer.data().clear();

    if (!showingUV) {
        viewer.data().set_mesh(V, F);
        viewer.data().set_face_based(false);

        if(UV.size() != 0) {
            viewer.data().set_uv(TextureResolution*UV);
            viewer.data().show_texture = true;
        }
    }
    else {
        viewer.data().show_texture = false;
        viewer.data().set_mesh(UV, F);
    }
	if (distortionMap.any())
		viewer.data().set_colors(distortionMap);
}

bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
	return false;
}


static void computeSurfaceGradientMatrix(SparseMatrix<double>& D1, SparseMatrix<double>& D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal() * Dx + F1.col(1).asDiagonal() * Dy + F1.col(2).asDiagonal() * Dz;
	D2 = F2.col(0).asDiagonal() * Dx + F2.col(1).asDiagonal() * Dy + F2.col(2).asDiagonal() * Dz;
}

void calculate_distortion() {

	Eigen::SparseMatrix<double> Dx, Dy;
	Eigen::MatrixXd J1, J2, J3, J4, J, D;

	distortionMap.conservativeResize(F.rows(), 3);
	J.conservativeResize(2, 2);
	distortion.conservativeResize(F.rows());

	computeSurfaceGradientMatrix(Dx, Dy);

	J1 = Dx * UV.col(0); J2 = Dx * UV.col(1); J3 = Dy * UV.col(0); J4 = Dy * UV.col(1);


	if (conformal)
	{
		Eigen::MatrixXd I = Eigen::MatrixXd::Identity(2, 2);

		for (int i = 0; i < F.rows(); i++) {
			J(0, 0) = J1(i, 0); J(0, 1) = J2(i, 0); J(1, 0) = J3(i, 0); J(1, 1) = J4(i, 0);

			D = J + J.transpose() - J.trace() * I;
			distortion[i] = pow(D.norm(), 2);
		}
	} else if (isometric) {
		Eigen::MatrixXd U, V, R, UVT;
		R.conservativeResize(2, 2);

		for (int i = 0; i < F.rows(); i++) {
			// Compute SVD
			J(0, 0) = J1(i, 0); J(0, 1) = J2(i, 0); J(1, 0) = J3(i, 0); J(1, 1) = J4(i, 0);
			JacobiSVD<MatrixXd> svd(J, ComputeThinU | ComputeThinV);
			U = svd.matrixU(); V = svd.matrixV(); UVT = U * V.transpose();
			R(0, 0) = 1; R(0, 1) = 0; R(1, 0) = 0; R(1, 1) = UVT.norm() == 0 ? 0 : 1;

			R = U * R * V.transpose();
			D = J - R;
			distortion[i] = pow(D.norm(), 2);
		}
	}
}

void show_distortion() {
	calculate_distortion();

	distortion = distortion / distortion.maxCoeff();

	for (int i = 0; i < F.rows(); i++) {
		distortionMap(i, 0) = 255;
		distortionMap(i, 1) = 1 - distortion[i];
		distortionMap(i, 2) = 1 - distortion[i];
	}

}

static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3)) * 0.5;
	double f = (J(0) - J(3)) * 0.5;
	double g = (J(1) + J(2)) * 0.5;
	double h = (J(1) - J(2)) * 0.5;
	double q = sqrt((e * e) + (h * h));
	double r = sqrt((f * f) + (g * g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1) * 0.5;
	double phi = (a2 + a1) * 0.5;

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

void ConvertConstraintsToMatrixForm(VectorXi indices, MatrixXd positions, Eigen::SparseMatrix<double>& C, VectorXd& d)
{
	// Convert the list of fixed indices and their fixed positions to a linear system
	// Hint: The matrix C should contain only one non-zero element per row and d should contain the positions in the correct order.
	Eigen::SparseMatrix<double> c;
	c.resize(indices.rows(), V.rows());
	C.resize(2 * indices.rows(), 2 * V.rows());

	for (int i = 0; i < indices.rows(); i++) {
		c.insert(i, indices(i, 0)) = 1;
	}

	igl::repdiag(c, 2, C);

	d.resize(2 * indices.rows(), 1);
	d = Map<VectorXd>(positions.data(), positions.cols() * positions.rows());
}

void computeParameterization(int type)
{
	VectorXi fixed_UV_indices;
	MatrixXd fixed_UV_positions;
	SparseMatrix<double> A, C;
	VectorXd b, d;

	// Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices
	if (!freeBoundary)
	{
		// The boundary vertices should be fixed to positions on the unit disc. Find these position and
		// save them in the #V x 2 matrix fixed_UV_position.
		igl::boundary_loop(F, fixed_UV_indices);
		igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
	}
	else
	{
		// Fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.
		vector<vector<int>> VV;
		igl::adjacency_list(F, VV);

		double max_distance = numeric_limits<double>::min();
		pair<int, int> indices;

		indices.first = -1; indices.second = -1;

		for (int i = 0; i < V.rows(); i++) {

			Eigen::VectorXd min_distance;
			Eigen::VectorXi previous;

			igl::dijkstra(i, {}, VV, min_distance, previous);

			for (int j = 0; j < min_distance.rows(); j++) {
				if (min_distance(j, 0) > max_distance) {
					max_distance = min_distance(j, 0);
					indices.first = i;
					indices.second = j;
				}
			}
		}

		fixed_UV_indices.resize(2, 1);

		fixed_UV_indices(0, 0) = indices.first;
		fixed_UV_indices(1, 0) = indices.second;

		igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);

	}

	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);

	// Find the linear system for the parameterization (1- Tutte, 2- Harmonic, 3- LSCM, 4- ARAP)
	// and put it in the matrix A.
	// The dimensions of A should be 2#V x 2#V.
	b.resize(2 * V.rows(), 1);
	b.setZero();
	if (type == '1') {
		// Add your code for computing uniform Laplacian for Tutte parameterization
		// Hint: use the adjacency matrix of the mesh
		Eigen::SparseMatrix<double> adjacency_mat, adjdiag, U;
		igl::adjacency_matrix(F, adjacency_mat);
		Eigen::SparseVector<double> adjsum;
		
		igl::sum(adjacency_mat, 1, adjsum);
		igl::diag(adjsum, adjdiag);
		U = adjacency_mat - adjdiag;
		igl::repdiag(U, 2, A);
	}

	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~
		Eigen::SparseMatrix<double> L;
		igl::cotmatrix(V, F, L);
		igl::repdiag(L, 2, A);
	}

	if (type == '3') {
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!
		Eigen::SparseMatrix<double> Dx, Dy, Dxt, Dyt, Jr, Jc, J, J1, J2, J3, J4;
		Eigen::VectorXd T;

		computeSurfaceGradientMatrix(Dx, Dy);
		igl::doublearea(V, F, T); T /= 2;
		
		auto Diag_T = T.asDiagonal();

		J1 = Dx.transpose() * Diag_T * Dx + Dy.transpose() * Diag_T * Dy;
		J2 = Dy.transpose() * Diag_T * Dx - Dx.transpose() * Diag_T * Dy;
		J3 = Dx.transpose() * Diag_T * Dy - Dy.transpose() * Diag_T * Dx;
		J4 = Dx.transpose() * Diag_T * Dx + Dy.transpose() * Diag_T * Dy;

		// Build J
		igl::cat(2, J1, J2, Jr);
		igl::cat(2, J3, J4, Jc);
		igl::cat(1, Jr, Jc, J);

		A = J;
	}

	// Solve the linear system.
	// Construct the system as discussed in class and the assignment sheet
	// Use igl::cat to concatenate matrices
	// Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail

	Eigen::SparseMatrix<double> N, M, Ct, lhs, zeroMat(Ct.cols(), C.rows());
	VectorXd rhs, res;

	Ct = C.transpose();

	igl::cat(1, A, C, N);
	igl::cat(1, Ct, zeroMat, M);
	igl::cat(2, N, M, lhs);

	rhs.resize(b.rows() + d.rows(), 1);
	rhs << b, d;

	Eigen::SparseLU<SparseMatrix<double>> solver;
	lhs.makeCompressed();
	solver.analyzePattern(lhs);
	solver.factorize(lhs);
	res = solver.solve(rhs);

	UV.resize(V.rows(), 2);

	for (int i = 0; i < V.rows(); i++) {
		UV(i, 0) = res(i, 0);
		UV(i, 1) = res(i + V.rows(), 0);
	}

}

bool callback_key_pressed(Viewer& viewer, unsigned char key, int modifiers) {
	switch (key) {
	case '1':
	case '2':
	case '3':
	case '4':
		computeParameterization(key);
		break;
	case '5':
	{
		show_distortion();
		break;
	}
	break;
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
		if (showingUV)
		{
			temp2D = viewer.core;
			viewer.core = temp3D;
			showingUV = false;
		}
		else
		{
			if (UV.rows() > 0)
			{
				temp3D = viewer.core;
				viewer.core = temp2D;
				showingUV = true;
			}
			else { std::cout << "ERROR ! No valid parameterization\n"; }
		}
		break;
	}
	Redraw();
	return true;
}

void init_core_states()
{
	// save initial viewer core state
	temp3D = viewer.core;
	temp2D = viewer.core;
	temp2D.orthographic = true;
}

bool load_mesh(string filename)
{
	igl::read_triangle_mesh(filename, V, F);
	Redraw();
	viewer.core.align_camera_center(V);
	showingUV = false;

	return true;
}

bool callback_init(Viewer& viewer)
{
	temp3D = viewer.core;
	temp2D = viewer.core;
	temp2D.orthographic = true;

	return false;
}

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cout << "Usage ex3_bin <mesh.off/obj>" << endl;
		load_mesh("../../../data/cathead.obj");
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
			static int chosen = 0;
			ImGui::RadioButton("Conformal Distortion Viszualization", &chosen, 0);
			ImGui::RadioButton("Isometric Distortion Viszualization", &chosen, 1);

			switch (chosen) {
			case 0:
				conformal = true;
				isometric = false;
				break;
			case 1:
				conformal = false;
				isometric = true;
				break;
			}

		}

		if (ImGui::Button("Upload Mesh", ImVec2(-1, 0)))
		{
			std::string fname = igl::file_dialog_open();
			std::cout << fname << " cloud has been uploaded\n";
			load_mesh(fname);
		}

		if (ImGui::Button("Init Core States", ImVec2(-1, 0)))
		{
			init_core_states();
		}
		
	};

	viewer.callback_key_pressed = callback_key_pressed;
	viewer.callback_mouse_move = callback_mouse_move;
	viewer.callback_init = callback_init;

	viewer.launch();
}
