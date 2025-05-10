#include "GCore/Components/MeshOperand.h"
#include "GCore/util_openmesh_bind.h"
#include "geom_node_base.h"
#include <cmath>
#include <Eigen/Sparse>

#include <pxr/base/gf/vec3f.h>
#include <pxr/base/vt/array.h>

/*
** @brief HW5_Laplacian_Surface_Editing
**
** This file presents the basic framework of a "node", which processes inputs
** received from the left and outputs specific variables for downstream nodes to
** use.
**
** - In the first function, node_declare, you can set up the node's input and
** output variables.
**
** - The second function, node_exec is the execution part of the node, where we
** need to implement the node's functionality.
**
**
** Your task is to fill in the required logic at the specified locations
** within this template, especially in node_exec.
*/

NODE_DEF_OPEN_SCOPE
NODE_DECLARATION_FUNCTION(mesh_editing)
{
    // Input-1: Original 3D mesh with boundary
	// Input-2: Position of all points after change
	// Input-3: Indices of control points

    b.add_input<Geometry>("Original mesh");
    b.add_input<pxr::VtVec3fArray>("Changed vertices");
    b.add_input<std::vector<size_t>>("Control Points Indices");

    /*
	** NOTE: You can add more inputs or outputs if necessary. For example, in
	** some cases, additional information (e.g. other mesh geometry, other
	** parameters) is required to perform the computation.
	**
	** Be sure that the input/outputs do not share the same name. You can add
	 one geometry as
	**
	**                b.add_input<Geometry>("Input");
	**
	*/


	// Output: New positions of all points or changed mesh
    b.add_output<pxr::VtVec3fArray>("New vertices");
}

NODE_EXECUTION_FUNCTION(mesh_editing)
{
	// Get the input from params
	auto input = params.get_input<Geometry>("Original mesh");
	pxr::VtVec3fArray changed_vertices = params.get_input<pxr::VtVec3fArray>("Changed vertices");
	std::vector<size_t> control_points_indices = params.get_input<std::vector<size_t>>("Control Points Indices");

	// Avoid processing the node when there is no input
	if (!input.get_component<MeshComponent>()) {
		throw std::runtime_error("Mesh Editing: Need Geometry Input.");
		return false;
	}
	if (control_points_indices.empty()) {
		throw std::runtime_error("at least one point");
	}

	/* ------------- [HW5_TODO] Laplacian Surface Editing Implementation -----------
	** Implement Laplacian surface editing to change the mesh.
	**
	** Steps:
	** 1. Building a large sparse matrix to solve the linear system.
	**
	** 2. Storage the solver and solve the linear system.The solver is used for pretreatment.
	**
	** 3. Update the positions of all points.
	**
	*/

	// Output the new positions of all points or you can output the changed mesh
	pxr::VtVec3fArray new_positions;
	params.set_output("New vertices", std::move(new_positions));

    return true;
}

NODE_DECLARATION_UI(mesh_editing);
NODE_DEF_CLOSE_SCOPE