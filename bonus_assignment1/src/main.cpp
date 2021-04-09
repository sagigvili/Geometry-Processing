
#include <igl/avg_edge_length.h>
#include <igl/read_triangle_mesh.h>
#include <igl/parula.h>
#include <igl/doublearea.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>


// Given (squared) edge-lengths of a triangle mesh `l_sqr` compute the internal
// angles at each corner (a.k.a. wedge) of the mesh.
//
// Inputs:
//    l_sqr  #F by 3 list squared edge-lengths opposite respective corner
// Outputs:
//    A  #F by 3 list of internal angles incident on respective corner
//
void internal_angles(
	const Eigen::MatrixXd& l_sqr,
	Eigen::MatrixXd& A)
{
	A.resizeLike(l_sqr);
	A.setConstant(1);
}

// Compute the discrete angle defect at each vertex of a triangle mesh
// (`V`,`F`), that is, the _locally integrated_ discrete Gaussian
// curvature.
//
// Inputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of mesh face indices into V
// Outputs:
//   D  #V list of angle defects in units radians
void angle_defect(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	Eigen::VectorXd& D)
{
	D.resize(V.rows());
	D.setConstant(1);
}

// Compute the discrete mean curvature at each vertex of a mesh (`V`,`F`) by
// taking the signed magnitude of the mean curvature normal as a _pointwise_ or
// _integral average_ quantity.
//
// Inputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of mesh face indices into V
// Outputs:
//   H  #V list of mean curvature values in units 1/m²
//
void mean_curvature(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	Eigen::VectorXd& H)
{
	H.resize(V.rows());
	H.setConstant(1);
}

// Approximate principal curvature values and directions locally by considering
// the two-ring neighborhood of each vertex in the mesh (`V`,`F`).
//
// Inputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of mesh face indices into V
// Outputs:
//   D1  #V by 3 list of first principal (unit) directions 
//   D2  #V by 3 list of second principal (unit) directions 
//   K1  #V by 3 list of first principal values (maximum curvature)
//   K2  #V by 3 list of second principal values (minimum curvature)
//
void principal_curvatures(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	Eigen::MatrixXd& D1,
	Eigen::MatrixXd& D2,
	Eigen::VectorXd& K1,
	Eigen::VectorXd& K2)
{
	D1.resizeLike(V);
	D2.resizeLike(V);
	K1.resize(V.rows());
	K2.resize(V.rows());
}

int main(int argc, char* argv[])
{
	// Scale for the color axis
	double scale = 100.0;
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	// Load input meshes
	igl::read_triangle_mesh(
		(argc > 1 ? argv[1] : "../data/cactus.obj"), V, F);

	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
	Eigen::VectorXd A = M.diagonal();

	Eigen::VectorXd D, G, H, K1, K2;
	Eigen::MatrixXd D1, D2;
	// Angle defect ~ locally integrated Gaussian curvature
	angle_defect(V, F, D);
	// average locally (i.e., "un-integrate" to pointwise quantity for
	// visualization)
	G = D.array() / A.array();
	mean_curvature(V, F, H);
	principal_curvatures(V, F, D1, D2, K1, K2);

	igl::opengl::glfw::Viewer viewer;
	std::cout << R"(
S,s      Stretch, squish color axis range
G        Show Gaussian curvature (using principal_curvatures)
g        Show Gaussian curvature (using angle_defect)
M        Show discrete mean curvature (using principal_curvatures)
m        Show discrete mean curvature (using mean_curvature)
K        Show maximum curvature (using principal_curvatures)
k        Show minimum curvature (using principal_curvatures)
D,d      Show principal directions 
)";
	// Default to mean curvature
	Eigen::VectorXd Z = H;
	const auto update = [&]()
	{
		Eigen::MatrixXd C;
		igl::parula(Z, -scale, scale, C);
		viewer.data().set_colors(C);
	};
	viewer.callback_key_pressed =
		[&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
	{
		switch (key)
		{
		case 'D':
		case 'd':
			viewer.data().show_overlay ^= 1;
			break;
		case 'G':
			Z = K1.array() * K2.array();
			break;
		case 'g':
			Z = G;
			break;
		case 'K':
			Z = K1;
			break;
		case 'k':
			Z = K2;
			break;
		case 'M':
			Z = 0.5 * (K1 + K2);
			break;
		case 'm':
			Z = H;
			break;
		case 'S':
		case 's':
			scale *= key == 'S' ? 2.0 : 0.5;
			std::cout << "Color axis range: [" << -scale << "," << scale << "]" << std::endl;
			break;
		default:
			return false;
		}
		update();
		return true;
	};

	viewer.data().set_mesh(V, F);
	Eigen::MatrixXd lP(V.rows() * 4, 3);
	const double h = igl::avg_edge_length(V, F);
	lP << V - 0.5 * h * D1, V + 0.5 * h * D1, V - 0.5 * h * D2, V + 0.5 * h * D2;
	Eigen::MatrixXi lE(2 * V.rows(), 2);
	Eigen::MatrixXd lC(2 * V.rows(), 3);

	const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
	const Eigen::RowVector3d yellow(1.0, 0.9, 0.2);
	const Eigen::RowVector3d blue(0.2, 0.3, 0.8);
	for (int e = 0;e < V.rows();e++)
	{
		lE(e, 0) = e + 0 * V.rows();
		lE(e, 1) = e + 1 * V.rows();
		lE(V.rows() + e, 0) = e + 2 * V.rows();
		lE(V.rows() + e, 1) = e + 3 * V.rows();
		lC.row(e) = orange;
		lC.row(V.rows() + e) = blue;
	}
	viewer.data().set_edges(lP, lE, lC);

	update();
	viewer.data().show_lines = false;
	viewer.data().show_overlay = false;
	viewer.data().face_based = false;
	viewer.launch();
	return EXIT_SUCCESS;
}
