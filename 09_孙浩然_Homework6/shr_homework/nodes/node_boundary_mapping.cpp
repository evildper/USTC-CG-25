#include "GCore/Components/MeshOperand.h"
#include "geom_node_base.h"
#include "GCore/util_openmesh_bind.h"
#include <Eigen/Sparse>

    /*
    ** @brief HW4_TutteParameterization
    **
    ** This file contains two nodes whose primary function is to map the boundary of
    *a mesh to a plain
    ** convex closed curve (circle of square), setting the stage for subsequent
    *Laplacian equation
    ** solution and mesh parameterization tasks.
    **
    ** Key to this node's implementation is the adept manipulation of half-edge data
    *structures
    ** to identify and modify the boundary of the mesh.
    **
    ** Task Overview:
    ** - The two execution functions (node_map_boundary_to_square_exec,
    ** node_map_boundary_to_circle_exec) require an update to accurately map the
    *mesh boundary to a and
    ** circles. This entails identifying the boundary edges, evenly distributing
    *boundary vertices along
    ** the square's perimeter, and ensuring the internal vertices' positions remain
    *unchanged.
    ** - A focus on half-edge data structures to efficiently traverse and modify
    *mesh boundaries.
    */

NODE_DEF_OPEN_SCOPE

    /*
    ** HW4_TODO: Node to map the mesh boundary to a circle.
    */

NODE_DECLARATION_FUNCTION(circle_boundary_mapping)
{
    // Input-1: Original 3D mesh with boundary
    b.add_input<Geometry>("Input");
    // Output-1: Processed 3D mesh whose boundary is mapped to a square and the
    // interior vertices remains the same
    b.add_output<Geometry>("Output");
}

NODE_EXECUTION_FUNCTION(circle_boundary_mapping)
{
    // Get the input from params
    auto input = params.get_input<Geometry>("Input");

    // (TO BE UPDATED) Avoid processing the node when there is no input
    // 
    if (!input.get_component<MeshComponent>()) {
        std::cerr << "Boundary Mapping: Need Geometry Input." << std::endl;
        return false;
    }
    
    

    /* ----------------------------- Preprocess -------------------------------
    ** Create a halfedge structure (using OpenMesh) for the input mesh. The
    ** half-edge data structure is a widely used data structure in geometric
    ** processing, offering convenient operations for traversing and modifying
    ** mesh elements.
    */
    auto halfedge_mesh = operand_to_openmesh(&input);
    //std::cout<<"sos"<<std::endl;
    /* ----------- [HW4_TODO] TASK 2.1: Boundary Mapping (to circle)
    *------------
    ** In this task, you are required to map the boundary of the mesh to a
    *circle
    ** shape while ensuring the internal vertices remain unaffected. This step
    *is
    ** crucial for setting up the mesh for subsequent parameterization tasks.
    **
    ** Algorithm Pseudocode for Boundary Mapping to Circle
    ** ------------------------------------------------------------------------
    ** 1. Identify the boundary loop(s) of the mesh using the half-edge
    *structure.
    **
    ** 2. Calculate the total length of the boundary loop to determine the
    *spacing
    **    between vertices when mapped to a square.
    **
    ** 3. Sequentially assign each boundary vertex a new position along the
    *square's
    **    perimeter, maintaining the calculated spacing to ensure proper
    *distribution.
    **
    ** 4. Keep the interior vertices' positions unchanged during this process.
    **
    ** Note: How to distribute the points on the circle?
    **
    ** Note: It would be better to normalize the boundary to a unit circle in
    *[0,1]x[0,1] for
    ** texture mapping.
    */

    /* ----------------------------- Postprocess ------------------------------
    ** Convert the result mesh from the halfedge structure back to Geometry
    *format as the node's
    ** output.
    */
   std::vector<OpenMesh::VertexHandle> boundary_vertices;
   OpenMesh::HalfedgeHandle he_start;
   
   // 1. 找到一个边界半边
   for (auto heh : halfedge_mesh->halfedges()) {
       if (halfedge_mesh->is_boundary(heh)) {
           he_start = heh;
           break;
       }
   }
   
   // 确保 he_start 是合法的
   if (!halfedge_mesh->is_boundary(he_start)) return false;
   
   // 2. 遍历边界顶点（按顺序存入 boundary_vertices）
   OpenMesh::HalfedgeHandle heh = he_start;
   do {
       boundary_vertices.push_back(halfedge_mesh->to_vertex_handle(heh));
       heh = halfedge_mesh->next_halfedge_handle(heh);  // **确保是顺时针遍历**
   } while (heh != he_start);
   
   // 3. 计算边界总长度
   double boundary_length = 0.0;
   for (size_t i = 0; i < boundary_vertices.size(); ++i) {
       OpenMesh::VertexHandle vh1 = boundary_vertices[i];
       OpenMesh::VertexHandle vh2 = boundary_vertices[(i + 1) % boundary_vertices.size()];
       boundary_length += (halfedge_mesh->point(vh1) - halfedge_mesh->point(vh2)).norm();
   }
   
   // 4. 映射到单位圆
   std::vector<OpenMesh::Vec3f> new_positions(boundary_vertices.size());
   double accumulated_length = 0.0;
   
   for (size_t i = 0; i < boundary_vertices.size(); ++i) {
       OpenMesh::VertexHandle vh = boundary_vertices[i];
   
       double theta = (accumulated_length / boundary_length) * 2.0 * M_PI;
       OpenMesh::Vec3d new_pos(std::cos(theta), std::sin(theta), 0.0);
       new_pos = OpenMesh::Vec3d(new_pos[0] * 0.5 + 0.5, new_pos[1] * 0.5 + 0.5, 0.0);
   
       new_positions[i] = OpenMesh::Vec3f(new_pos[0], new_pos[1], new_pos[2]);
   
       OpenMesh::VertexHandle next_vh = boundary_vertices[(i + 1) % boundary_vertices.size()];
       accumulated_length += (halfedge_mesh->point(vh) - halfedge_mesh->point(next_vh)).norm();
   }
   
   // **5. 统一更新顶点位置**
   for (size_t i = 0; i < boundary_vertices.size(); ++i) {
       halfedge_mesh->set_point(boundary_vertices[i], new_positions[i]);
   }
   
   // 6. 输出处理后的网格
   auto geometry = openmesh_to_operand(halfedge_mesh.get());
   params.set_output("Output", std::move(*geometry));
   
   return true;
   
}

    /*
    ** HW4_TODO: Node to map the mesh boundary to a square.
    */

NODE_DECLARATION_FUNCTION(square_boundary_mapping)
{
    // Input-1: Original 3D mesh with boundary
    b.add_input<Geometry>("Input");

    // Output-1: Processed 3D mesh whose boundary is mapped to a square and the
    // interior vertices remains the same
    b.add_output<Geometry>("Output");
}

NODE_EXECUTION_FUNCTION(square_boundary_mapping)
{
    // Get the input from params
    auto input = params.get_input<Geometry>("Input");

    // (TO BE UPDATED) Avoid processing the node when there is no input
    // if (!input.get_component<MeshComponent>()) {
    //     throw std::runtime_error("Input does not contain a mesh");
    // }
    // throw std::runtime_error("Not implemented");
    if (!input.get_component<MeshComponent>()) {
        std::cerr << "Boundary Mapping: Need Geometry Input." << std::endl;
        return false;
    }
    /* ----------------------------- Preprocess -------------------------------
    ** Create a halfedge structure (using OpenMesh) for the input mesh.
    */
    auto halfedge_mesh = operand_to_openmesh(&input);

    /* ----------- [HW4_TODO] TASK 2.2: Boundary Mapping (to square)
    *------------
    ** In this task, you are required to map the boundary of the mesh to a
    *circle
    ** shape while ensuring the internal vertices remain unaffected.
    **
    ** Algorithm Pseudocode for Boundary Mapping to Square
    ** ------------------------------------------------------------------------
    ** (omitted)
    **
    ** Note: Can you perserve the 4 corners of the square after boundary
    *mapping?
    **
    ** Note: It would be better to normalize the boundary to a unit circle in
    *[0,1]x[0,1] for
    ** texture mapping.
    */

    /* ----------------------------- Postprocess ------------------------------
    ** Convert the result mesh from the halfedge structure back to Geometry
    *format as the node's
    ** output.
    */
    std::vector<OpenMesh::VertexHandle> boundary_vertices;
    OpenMesh::HalfedgeHandle he_start;
  
  // 找到一个边界半边
    for (auto heh : halfedge_mesh->halfedges()) {
        if (halfedge_mesh->is_boundary(heh)) {
            he_start = heh;
            break;
        }
    }
    // 确保 he_start 有效
    if (!he_start.is_valid()) return false;

  // 沿着边界遍历，按顺序存入 boundary_vertices
    OpenMesh::HalfedgeHandle heh = he_start;
    do {
        boundary_vertices.push_back(halfedge_mesh->to_vertex_handle(heh));
        heh = halfedge_mesh->next_halfedge_handle(heh);  // 沿着边界移动到下一个半边
    } while (heh != he_start);
    double total_length = 0.0;
    for (size_t i = 0; i < boundary_vertices.size(); ++i) {
        OpenMesh::VertexHandle vh = boundary_vertices[i];
        OpenMesh::VertexHandle next_vh = boundary_vertices[(i + 1) % boundary_vertices.size()];
        total_length += (halfedge_mesh->point(vh) - halfedge_mesh->point(next_vh)).norm();
    }
    double accumulated_length = 0.0;//通过边长比来映射到[0,1]×[0,1]正方形边界
    for (size_t i =0;i<boundary_vertices.size();i++)
    {   
        double t = accumulated_length/total_length;
        OpenMesh::VertexHandle vh = boundary_vertices[i];
        OpenMesh::VertexHandle next_vh = boundary_vertices[(i + 1) % boundary_vertices.size()];
        accumulated_length += (halfedge_mesh->point(vh) - halfedge_mesh->point(next_vh)).norm();
        OpenMesh::Vec3d new_pos(0.0, 0.0, 0.0);
        if(t<0.25)
        {
            new_pos = OpenMesh::Vec3d(4.0*t,0.0,0.0);
        }
        else if(t<0.5)
        {
            new_pos = OpenMesh::Vec3d(1.0,4.0*(t-0.25),0.0);
        }
        else if(t<0.75)
        {
            new_pos = OpenMesh::Vec3d(1.0-4.0*(t-0.5),1.0,0.0);
        }
        else
        {
            new_pos = OpenMesh::Vec3d(0.0,1.0-4.0*(t-0.75),0.0);
        }
        halfedge_mesh->set_point(vh, OpenMesh::Vec3f(new_pos[0], new_pos[1], new_pos[2]));
    }
    auto geometry = openmesh_to_operand(halfedge_mesh.get());

    // Set the output of the nodes
    params.set_output("Output", std::move(*geometry));
    return true;
}



NODE_DECLARATION_UI(boundary_mapping);
NODE_DEF_CLOSE_SCOPE