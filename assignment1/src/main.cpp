#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/writeOFF.h>

#include "gui_utils.h"

/*** insert any libigl headers here ***/
#include <igl/facet_components.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/boundary_loop.h>
#include <igl/loop.h>
#include <igl/edge_topology.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>

using namespace std;

enum MouseMode { NONE, FACE_SELECT, VERTEX_SELECT, TRANSLATE, ROTATE };
MouseMode mouse_mode = NONE;
bool doit = false;
int down_mouse_x = -1, down_mouse_y = -1;
// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd colors_per_face;

std::set<int> selected_faces;
int selected_v = -1;
void update_display(igl::opengl::glfw::Viewer& viewer);

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
    if (key == '1')
    {
        cout << endl;
        igl::vertex_triangle_adjacency(V, F, VF, VFi);
        cout << "Vertex        Adjacented Faces";
        cout << endl;
        for (auto VFit = VF.begin(); VFit != VF.end(); VFit++)
        {
            cout << "-------------------------------" << endl;
            cout << "  " << VFit - VF.begin() << "            ";
            std::vector<int> adjs = *VFit;
            for (auto AJit = adjs.begin(); AJit != adjs.end(); AJit++)
            {
                cout << *AJit << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

	if (key == '2')
	{
        cout << endl;
        igl::adjacency_list(F, VV);
        cout << "Vertex        Adjacented Vertices";
        cout << endl;
        for (auto VVit = VV.begin(); VVit != VV.end(); VVit++)
        {
            cout << "-------------------------------" << endl;
            cout << "  " << VVit - VV.begin() << "            ";
            std::vector<int> adjs = *VVit;
            for (auto AJit = adjs.begin(); AJit != adjs.end(); AJit++)
            {
                cout << *AJit << " ";
            }
            cout << endl;
        }
        cout << endl;
	}

	if (key == '3')
	{
		viewer.data().clear();
		viewer.data().set_mesh(V, F);
		colors_per_face.setZero(F.rows(),3);
        igl::facet_components(F, cid);
        igl::jet(cid, true, colors_per_face);

        // cid contains now every face and the component number it's part of
        // maxCoeff gets the highest number in a VectorXi, cid in this case
        // The componenets starts from 0, so we add 1
        int num_of_comps = cid.maxCoeff() + 1;
        cout << endl << "Number of componenets: " << num_of_comps << endl;
        cout << "Component        Number of faces" << endl;
        
        for (int i = 0; i < num_of_comps; i++) {
            int faces = 0;
            for (auto CIDit = cid.begin(); CIDit != cid.end(); CIDit++)
            {
                if (i == *CIDit)
                    faces++;
            }
            cout << i << "                " << faces << endl;
        }
		viewer.data().set_colors(colors_per_face);
	}

	if (key == '4')
	{
        double a, n;
        int count = 0;
		Eigen::MatrixXd Vout=V , P, baryCent;
		Eigen::MatrixXi Fout=F, M;
        Fout.setZero(F.rows() * 3, F.cols());
        P.setZero(V.rows(), V.cols());

        igl::triangle_triangle_adjacency(F, M);
        igl::barycenter(V, F, baryCent);
        igl::adjacency_list(F, VV);

        for (int i = 0; i < V.rows(); i++)
        {
            n = VV[i].size();
            float piX2 = 3.14159265 * 2;
            a = (4 - 2 * cos(piX2 / n)) / 9;
            for (int j = 0; j < VV[i].size(); j++)
            {
                P.row(i) += V.row(VV[i][j]);
            }
            P.row(i) = ((1 - a) * V.row(i)) + (a * (P.row(i) / n));
        }
        Vout.setZero(P.rows() + baryCent.rows(), P.cols());
        Vout << P, baryCent;

        for (int i = 0; i < F.rows(); i++)
        {
            for (int j = 0; j < 3; j++) {
                if (M(i, j) == -1)
                    Fout.row(count++) << F(i, j), F(i, (j + 1) % 3), V.rows() + i;
                else if (i > M(i, j))
                {
                    Fout.row(count) << F(i, j), V.rows() + M(i, j), V.rows() + i;
                    Fout.row(++count) << V.rows() + M(i, j), F(i, (j + 1) % 3), V.rows() + i;
                    count++;
                }
            }
        }

		V = Vout; 
        F = Fout;
		update_display(viewer);
	}
    
    return false;
}

std::set<int> get_v_from_faces_idx(const Eigen::MatrixXi& F, std::set<int>& face_idx) {
    std::set<int> v_set;
    for (auto f: face_idx) {
        v_set.insert(F(f,0)); v_set.insert(F(f,1)); v_set.insert(F(f,2));
    }
    return v_set;
}

void extrude(igl::opengl::glfw::Viewer& viewer) {
    Eigen::MatrixXd Vout=V;
    Eigen::MatrixXi Fout=F;

    // Get selected faces
    Eigen::MatrixXi sF(selected_faces.size(),3); int idx = 0;
    for (auto it=selected_faces.begin();it!=selected_faces.end();it++){sF.row(idx++) = F.row(*it);}

    // Assert selected faces are connected
    Eigen::VectorXi comp; igl::facet_components(sF,comp);
    if (comp.maxCoeff() != 0) { cout << "Error: Not a single connected component, #face_comp =  " << comp << endl; return;}

    // 1) Get the boundary vertices surrounding the selected faces
    std::vector<int> bnd_loop; 
    igl::boundary_loop(sF,bnd_loop);

    // 2) Duplicate boundary vertices
    Vout.resize(V.rows()+bnd_loop.size(),3);
    for (int i = 0; i < V.rows(); i++) Vout.row(i)=V.row(i); // set vertices as old vertices
    for (int i = 0; i < bnd_loop.size(); i++) {Vout.row(V.rows()+i) = V.row(bnd_loop[i]);} // create new vertices as duplicates of the faces boundary

    // 3) Compute direction T: The average of selected face normals
    Eigen::RowVector3d T; T.setZero(); Eigen::MatrixXd FN;
    igl::per_face_normals(V,F,FN);
    for (auto it=selected_faces.begin();it!=selected_faces.end();it++){ T += FN.row(*it);}
    T.normalize(); T*=0.25*(V.row(bnd_loop[1])-V.row(bnd_loop[0])).norm();

    // 4) Offset old vertices by T
    std::set<int> inner_v = get_v_from_faces_idx(F, selected_faces);;
    for (auto v: inner_v) {
        Vout.row(v) += T;
    }

    // 5) Update Fout 
    Fout.resize(F.rows()+2*bnd_loop.size(),3); // 2 new faces per new edge (= per new vertex)
    for (int i = 0; i < F.rows(); i++) Fout.row(i)=F.row(i); // set first 'F.rows()' faces as the old faces


    // 5.1) Get the set of faces containing the old boundary vertices (hint: call igl::vertex_triangle_adjacency on the old 'F')
    igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);
    
    // 5.2) Get the "outer" set of faces containing the boundary vertices 
    //      (hint: call std::set_difference to compute the difference between the previously computed set of faces, and the selected faces)
    std::set<int> boundery_verices, outer_faces;
    for (auto BNDit = bnd_loop.begin(); BNDit != bnd_loop.end(); BNDit++)
    {
        int vf_size = VF[*BNDit].size();
        for (int j = 0; j < vf_size; j++) {
            boundery_verices.insert(VF[*BNDit][j]);
        }
    }
    std::set_difference(boundery_verices.begin(), boundery_verices.end(), selected_faces.begin(), selected_faces.end(), inserter(outer_faces, outer_faces.end()));

    // 5.3) Edit old outer faces indices, replacing the old vertices with the indices of the duplicated boundary vertices
    for (auto OUTERit = outer_faces.begin(); OUTERit != outer_faces.end(); OUTERit++)
        for (int i = 0; i < bnd_loop.size(); i++)
            for (int j = 0; j < 3; j++)
                if (Fout.row(*OUTERit)[j] == bnd_loop[i])
                    Fout.row(*OUTERit)[j] = V.rows() + i;

    // 5.4) Add new faces, 2 per edge
    int f_idx = F.rows();
    for (int i = 0; i < bnd_loop.size(); i++) {
        int v1,v2,v3,v4, j = (i + 1) % bnd_loop.size();
        v1 = V.rows() + i;
        v2 = V.rows() + j;
        v3 = bnd_loop[j];
        v4 = bnd_loop[i];
        Fout.row(f_idx++) << v1,v2,v3;
        Fout.row(f_idx++) << v3,v4,v1;
    }

    // 6) Check that the new mesh is a manifold (call is_edge_manifold, is_vertex_manifold on Vout,Fout)
    Eigen::VectorXi t1, t2;
    igl::is_vertex_manifold(Fout, t1);
    t2.setOnes(t1.rows(), t1.cols());
    if (!igl::is_edge_manifold(Fout) || t1 != t2) 
        return;

    // 7) Update V,F
    V = Vout;
    F = Fout;

    // Update gui and move to edit-translate mode
    colors_per_face = Eigen::MatrixXd::Ones(F.rows(),3); // number of faces has changed
    viewer.data().clear();
    viewer.data().set_mesh(V,F);
    for (auto f: selected_faces) {colors_per_face.row(f) << 1,0,0;}
    viewer.data().set_colors(colors_per_face);
    mouse_mode = TRANSLATE;
}

void clear_selection(igl::opengl::glfw::Viewer& viewer) {
    selected_faces.clear();
    selected_v = -1;
    colors_per_face = Eigen::MatrixXd::Ones(F.rows(),3);
    viewer.data().clear();
    viewer.data().set_mesh(V,F);
    viewer.data().set_colors(colors_per_face);
}

void export_mesh() {
    std::string f = igl::file_dialog_save();
    igl::writeOFF(f,V,F);
}

bool callback_init(igl::opengl::glfw::Viewer& viewer)
{
    return false;
}

bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) 
{
    down_mouse_x = viewer.current_mouse_x;
    down_mouse_y = viewer.current_mouse_y;

    if (mouse_mode == FACE_SELECT) 
	{
        int f = pick_face(viewer, down_mouse_x, down_mouse_y,V,F);
        if (f !=-1)  
		{
            selected_faces.insert(f);
            selected_v = -1;
            // update face colors
            //colors_per_face.setConstant(colors_per_face.rows(), colors_per_face.cols(), 0);
            colors_per_face.row(f) << 1,0,0;
            viewer.data().set_colors(colors_per_face);
        }
        
    } 
	else if (mouse_mode == VERTEX_SELECT) 
	{
        int v = pick_vertex(viewer, down_mouse_x, down_mouse_y,V,F);
        if (v !=-1) 
		{
            selected_v = v;
            selected_faces.clear(); 
            update_display(viewer);
            viewer.data().set_points(V.row(selected_v),Eigen::RowVector3d(1,0,0));
        }
    } 
	else if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE)) 
	{
        if (!selected_faces.empty()) 
		{
            int f = pick_face(viewer, down_mouse_x, down_mouse_y,V,F);
            if (std::find(selected_faces.begin(),selected_faces.end(),f)!= selected_faces.end()) 
			{
                doit = true;
            }
        } 
		else if (selected_v != -1) 
		{
            int v = pick_vertex(viewer, down_mouse_x, down_mouse_y,V,F);
            if (v == selected_v) 
			{
                viewer.data().set_points(V.row(selected_v),Eigen::RowVector3d(1,0,0));
                doit = true;
            }
        }
    }
    return false;
}

Eigen::RowVector3d get_face_avg(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, std::set<int>& selected_faces) {
    
    Eigen::RowVector3d avg; avg << 0,0,0;
    std::set<int> v_set = get_v_from_faces_idx(F, selected_faces);
    for (auto v: v_set) {
        avg += V.row(v);
    }
    avg/= v_set.size();
    
    return avg;
}

bool callback_mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y) 
{
	if (!doit)
	{
		return false;
	}
    if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
	{
        if (!selected_faces.empty())
		{
            Eigen::RowVector3d face_avg_pt = get_face_avg(V,F,selected_faces);
            std::set<int> v_idx = get_v_from_faces_idx(F,selected_faces);
            if (mouse_mode == TRANSLATE)
			{
                Eigen::Vector3f translation = computeTranslation(viewer,
                                             mouse_x,
                                             down_mouse_x,
                                             mouse_y,
                                             down_mouse_y,
                                             face_avg_pt);

                for (auto v_i : v_idx) 
				{
					V.row(v_i) += translation.cast<double>();
				}
            } 
			else 
			{ // ROTATE
                Eigen::Vector4f rotation = computeRotation(viewer,
                                 mouse_x,
                                 down_mouse_x,
                                 mouse_y,
                                 down_mouse_y,
                                 face_avg_pt);

                for (auto v_i : v_idx) 
				{
					Eigen::RowVector3f goalPosition = V.row(v_i).cast<float>();
					goalPosition -= face_avg_pt.cast<float>();
					Eigen::RowVector4f goalPosition4 = Eigen::RowVector4f(goalPosition[0], goalPosition[1], goalPosition[2], 0);
					igl::rotate_by_quat(goalPosition4.data(), rotation.data(), goalPosition4.data());
					goalPosition = Eigen::RowVector3f(goalPosition4[0], goalPosition4[1], goalPosition4[2]);
					goalPosition += face_avg_pt.cast<float>();
					V.row(v_i) = goalPosition.cast<double>();
                }
            }

            viewer.data().set_mesh(V,F);
            down_mouse_x = mouse_x;
            down_mouse_y = mouse_y;
            return true;    
        } 
		else if ((selected_v!=-1) && (mouse_mode == TRANSLATE)) 
		{
            Eigen::Vector3f translation = computeTranslation(viewer,
                                             mouse_x,
                                             down_mouse_x,
                                             mouse_y,
                                             down_mouse_y,
                                             V.row(selected_v));

            V.row(selected_v) += translation.cast<double>();
            viewer.data().set_mesh(V,F);
            viewer.data().set_points(V.row(selected_v),Eigen::RowVector3d(1,0,0));
            down_mouse_x = mouse_x;
            down_mouse_y = mouse_y;
            return true;
        }
    }

    return false;
}

bool callback_mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
    doit = false;
    return false;
}

void update_display(igl::opengl::glfw::Viewer& viewer) {
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    colors_per_face = Eigen::MatrixXd::Ones(F.rows(),3);
    viewer.data().set_colors(colors_per_face);
}

int main(int argc, char *argv[]) {
	// Init the viewer
	igl::opengl::glfw::Viewer viewer;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	// Add content to the default menu window
	menu.callback_draw_viewer_menu = [&]()
	{
		menu.draw_viewer_menu();

		ImGui::Combo("Mouse Mode", (int *)(&mouse_mode), "None\0Face Selection\0Vertex Selection\0Translate\0Roate\0\0");

		if (ImGui::Button("Extrude"))
		{
			extrude(viewer);
		}

		if (ImGui::Button("Clear Selection"))
		{
			clear_selection(viewer);
		}

		if (ImGui::Button("Export Mesh"))
		{
			export_mesh();
		}
	};

    //Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    
    if (argc == 2)
    {
      // Read mesh
      igl::readOFF(argv[1],V,F);
      
    }
    else
    {
      // Read mesh
      igl::readOFF("../data/cube.off",V,F);
    }

    viewer.data().set_mesh(V,F);
    viewer.data().compute_normals();
    viewer.core.align_camera_center(V,F);
    viewer.callback_init = callback_init;
    viewer.callback_mouse_down = callback_mouse_down;
    viewer.callback_mouse_up = callback_mouse_up;
    viewer.callback_mouse_move = &callback_mouse_move;
    update_display(viewer);
    viewer.launch();
}
