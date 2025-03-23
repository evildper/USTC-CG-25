#include "GCore/Components/MeshOperand.h"
#include "GCore/util_openmesh_bind.h"
#include "geom_node_base.h"
#include <cmath>
#include <time.h>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Geometry>

    /*
    ** @brief HW4_TutteParameterization
    **
    ** This file presents the basic framework of a "node", which processes inputs
    ** received from the left and outputs specific variables for downstream nodes to
    ** use.
    ** - In the first function, node_declare, you can set up the node's input and
    ** output variables.
    ** - The second function, node_exec is the execution part of the node, where we
    ** need to implement the node's functionality.
    ** Your task is to fill in the required logic at the specified locations
    ** within this template, especially in node_exec.
    */

NODE_DEF_OPEN_SCOPE
NODE_DECLARATION_FUNCTION(min_surf2)//写在node声明之下，实现之前。
{
    // Input-1: Original 3D mesh with boundary
    b.add_input<Geometry>("Input");

    /*
    ** NOTE: You can add more inputs or outputs if necessary. For example, in
    *some cases,
    ** additional information (e.g. other mesh geometry, other parameters) is
    *required to perform
    ** the computation.
    **
    ** Be sure that the input/outputs do not share the same name. You can add
    *one geometry as
    **
    **                b.add_input<Geometry>("Input");
    **
    ** Or maybe you need a value buffer like:
    **
    **                b.add_input<float1Buffer>("Weights");
    */

    // Output-1: Minimal surface with fixed boundary
    b.add_output<Geometry>("Output");
}
double compute_cotangent_weight(OpenMesh::HalfedgeHandle he, OpenMesh::PolyMesh_ArrayKernelT<>& mesh) {
    if (mesh.is_boundary(he)) return 0.0;

    auto vi = mesh.from_vertex_handle(he);  // 顶点 i
    auto vj = mesh.to_vertex_handle(he);    // 顶点 j
    auto opp_he = mesh.opposite_halfedge_handle(he);// 把边反向

    double alpha = 0.0, beta = 0.0;

    // 计算 α_ij（第一侧的对角点）
    if (!mesh.is_boundary(he)) {
        auto vk = mesh.to_vertex_handle(mesh.next_halfedge_handle(he));  // 找到对角点 vk , 同一个三角形半边的方向是一致的 。
        Eigen::Vector3d a = Eigen::Vector3d(mesh.point(vi)[0], mesh.point(vi)[1], mesh.point(vi)[2]) - Eigen::Vector3d(mesh.point(vk)[0], mesh.point(vk)[1], mesh.point(vk)[2]);
        //注意mesh.point(vi)是OpenMesh:vector 3f 类型 ， 无法直接赋值给Eigen::vector 3f类型。
        Eigen::Vector3d b = Eigen::Vector3d(mesh.point(vj)[0], mesh.point(vj)[1], mesh.point(vj)[2]) - Eigen::Vector3d(mesh.point(vk)[0], mesh.point(vk)[1], mesh.point(vk)[2]);

        double cross_norm = a.cross(b).norm();//做叉乘再求模长
        alpha = std::atan2(cross_norm,a.dot(b));//利用arctan求弧度
    }

    // 计算 β_ij（另一侧的对角点）
    if (!mesh.is_boundary(opp_he)) {
        auto vl = mesh.to_vertex_handle(mesh.next_halfedge_handle(opp_he));  // 找到对角点 vl
        Eigen::Vector3d a = Eigen::Vector3d(mesh.point(vi)[0], mesh.point(vi)[1], mesh.point(vi)[2]) - Eigen::Vector3d(mesh.point(vl)[0], mesh.point(vl)[1], mesh.point(vl)[2]);

        Eigen::Vector3d b = Eigen::Vector3d(mesh.point(vj)[0], mesh.point(vj)[1], mesh.point(vj)[2]) - Eigen::Vector3d(mesh.point(vl)[0], mesh.point(vl)[1], mesh.point(vl)[2]);

        double cross_norm = a.cross(b).norm();
        beta = std::atan2(cross_norm,a.dot(b));
    }
    double length = (mesh.point(vi)-mesh.point(vj)).norm();
    double weight = (std::tan(alpha/2.)+std::tan(beta/2.))/length;
    return std::max(weight,0.); //避免返回负的权重。                                                                                                                                                                                                                                                                                                                                                                                                  
}

NODE_EXECUTION_FUNCTION(min_surf2)
{
    auto input = params.get_input<Geometry>("Input");
    if (!input.get_component<MeshComponent>()) {
        throw std::runtime_error("Minimal Surface: Need Geometry Input.");
        return false;
    }
    auto halfedge_mesh = operand_to_openmesh(&input);//std::shared_ptr<USTC_CG::PolyMesh> halfedge_mesh
    size_t  num_vertices = halfedge_mesh->n_vertices();//计算顶点个数，确定矩阵大小
    Eigen::SparseMatrix<double> A(num_vertices, num_vertices);  // 拉普拉斯矩阵
    std::vector<Eigen::VectorXd> B_channels(3, Eigen::VectorXd( num_vertices ));
    std::vector<Eigen::Triplet<double>> triplet_list;
    for (auto vh : halfedge_mesh->vertices()) {
        

        if(halfedge_mesh->is_boundary(vh))
        {
            triplet_list.emplace_back(vh.idx(),vh.idx(),1.0);
            for(size_t i = 0; i<3; i++)
            {
                B_channels[i](vh.idx()) = halfedge_mesh->point(vh)[i];//半边数据结构本身对点有个idx(),从0，1，2.....
            }

        }
        else
        {
            double sum = 0;
            for(auto v : halfedge_mesh->vv_range(vh))//vv_range(vh):遍历周围的点
            {   
                auto he = halfedge_mesh->find_halfedge(vh, v);//找到半边
                double w = compute_cotangent_weight(he, *halfedge_mesh);//获得此半边对应的权重
                if(w > 0.)//避免负值插入
                {
                    triplet_list.emplace_back(vh.idx(),v.idx(),-w);
                    sum +=w;
                }
                
            }
            triplet_list.emplace_back(vh.idx(),vh.idx(),sum);//归一化
            for(size_t i = 0; i<3; i++)
            {
                B_channels[i](vh.idx()) = 0.0;
            }
        }

    }
    A.setFromTriplets(triplet_list.begin(), triplet_list.end());
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    std::vector<Eigen::VectorXd> X_channels(3);
    for (int c = 0; c < 3; ++c) {
        X_channels[c] = solver.solve(B_channels[c]);  
    }
    for (auto v : halfedge_mesh->vertices()) {
        if (!halfedge_mesh->is_boundary(v)) {
            halfedge_mesh->set_point(v, OpenMesh::Vec3f(X_channels[0](v.idx()),X_channels[1](v.idx()),X_channels[2](v.idx())));            
        }
    }   //这里也要转换成OpenMesh::Vec3格式
    auto geometry = openmesh_to_operand(halfedge_mesh.get());

    // Set the output of the nodes
    params.set_output("Output", std::move(*geometry));
    return true;
}

NODE_DECLARATION_UI(min_surf2);
NODE_DEF_CLOSE_SCOPE
