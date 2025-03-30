#include "GCore/Components/MeshOperand.h"
#include "GCore/util_openmesh_bind.h"
#include "geom_node_base.h"
#include <cmath>
#include <time.h>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <Eigen/Sparse>
#include "arap.h"

/*
** @brief HW5_ARAP_Parameterization
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
** - The third function generates the node's registration information, which
** eventually allows placing this node in the GUI interface.
**
** Your task is to fill in the required logic at the specified locations
** within this template, especially in node_exec.
*/

NODE_DEF_OPEN_SCOPE
NODE_DECLARATION_FUNCTION(arap)
{
    // Input-1: Original 3D mesh with boundary
    // Maybe you need to add another input for initialization?
    b.add_input<Geometry>("Input");
    b.add_input<Geometry>("Origin Mesh");
    
    b.add_input<int>("fix index1").min(0).max(3000).default_val(0);
    b.add_input<int>("fix index2").min(0).max(3000).default_val(0);


    /*
    ** NOTE: You can add more inputs or outputs if necessary. For example, in
    ** some cases, additional information (e.g. other mesh geometry, other
    ** parameters) is required to perform the computation.
    **
    ** Be sure that the input/outputs do not share the same name. You can add
    ** one geometry as
    **
    **                b.add_input<Geometry>("Input");
    **
    ** Or maybe you need a value buffer like:
    **
    **                b.add_input<float1Buffer>("Weights");
    */

    // Output-1: The UV coordinate of the mesh, provided by ARAP algorithm
    b.add_output<pxr::VtArray<pxr::GfVec2f>>("OutputUV");
    b.add_output<Geometry>("Output");
}

NODE_EXECUTION_FUNCTION(arap)
{
    // Get the input from params
    auto input = params.get_input<Geometry>("Input");
    auto input2 = params.get_input<Geometry>("Origin Mesh");
    // Avoid processing the node when there is no input
    // if (!(input.get_component<MeshComponent>() && input2.get_component<MeshComponent>())) {
    //     throw std::runtime_error("Need Geometry Input.");
    // }
    int fix1 = params.get_input<int>("fix index1");
    int fix2 = params.get_input<int>("fix index2");
    /* ----------------------------- Preprocess -------------------------------
    ** Create a halfedge structure (using OpenMesh) for the input mesh. The
    ** half-edge data structure is a widely used data structure in geometric
    ** processing, offering convenient operations for traversing and modifying
    ** mesh elements.
    */
    auto halfedge_mesh = operand_to_openmesh(&input);
    auto origin_mesh = operand_to_openmesh(&input2);

    Arap* arap = new Arap;
    double energy = 0, energy_new = 1;
    arap->init(origin_mesh, halfedge_mesh);
    arap->set_uv_mesh();
    arap->set_flatxy();
    arap->set_fixed(fix1, fix2);
    arap->set_cotangent();
    arap->set_matrixA();

    /*
    arap->set_matrixLt();
    arap->SVD_Lt();
    arap->set_Laplacian();
    energy_new = arap->energy_cal();
    arap->set_new_mesh();
    arap->reset_mesh();
    */
    
    int flag = 0;
    while (fabs(energy - energy_new) > 0.01 && flag < 100)
    { 
        energy = energy_new;
        arap->set_matrixLt();
        arap->SVD_Lt();
        arap->set_Laplacian();
        energy_new = arap->energy_cal();
        arap->set_new_mesh();
        arap->reset_mesh();
        flag++;
    }
    
    // The result UV coordinates
    pxr::VtArray<pxr::GfVec2f> uv_result;
    uv_result = arap->get_new_mesh();

    // Set the output of the node
    //
    for (const auto& vertex_handle : halfedge_mesh->vertices())
    {
        const auto& vec = halfedge_mesh->point(vertex_handle);
        halfedge_mesh->set_point(vertex_handle, OpenMesh::Vec3f(uv_result[vertex_handle.idx()][0],uv_result[vertex_handle.idx()][1] , 0));
        
    }
    //
    //params.set_output("OutputUV", uv_result);
    auto geometry = openmesh_to_operand(halfedge_mesh.get());
    
    // // Set the output of the nodes
    params.set_output("OutputUV", uv_result);
    params.set_output("Output", std::move(*geometry));
    
    return true;
}

NODE_DECLARATION_UI(arap);
NODE_DEF_CLOSE_SCOPE