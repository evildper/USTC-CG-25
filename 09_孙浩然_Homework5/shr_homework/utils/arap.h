#pragma once

#include "GCore/Components/MeshOperand.h"
#include "GCore/util_openmesh_bind.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace USTC_CG
{
	class Arap
	{
		public:
            Arap() = default;

			~Arap() = default;
            // 输入初变形和原始的半边结构
            void init(std::shared_ptr<PolyMesh> mesh1, std::shared_ptr<PolyMesh> mesh2);
			// 计算每条边的权重
			void set_cotangent();
			// 初始化矩阵
			void set_matrixA();
            // 设置两个固定点
            void set_fixed(int a, int b);
			// 得到初始uv坐标
            void set_uv_mesh();
			// 将空间的三角形投影到二维平面上
            void set_flatxy();
			// 计算每个面的仿射矩阵
			void set_matrixLt();
			// 进行SVD分解
            void SVD_Lt();
			// 设置拉普拉斯
            void set_Laplacian();
			// 计算能量函数
            double energy_cal();
			// 得到新的半边结构
            void set_new_mesh();
			// 存储半边结构
            void reset_mesh();
			// 得到新的uv坐标。
            pxr::VtArray<pxr::GfVec2f> get_new_mesh();

		private:
            //两个半边结构
            std::shared_ptr<PolyMesh> origin_mesh, mesh;
            //每条边的参数
            pxr::VtArray<double> cotangent;
            //稀疏矩阵A
            Eigen::SparseMatrix<double> A;
            //新旧坐标
            Eigen::VectorXd bx, by, new_x, new_y;
            //稀疏求解器
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            //每个面的仿射矩阵
            pxr::VtArray<Eigen::Matrix2d> Lt;
            //每个面的旋转矩阵
            pxr::VtArray<Eigen::Matrix2d> SVDLt;
            //新旧uv坐标
            pxr::VtArray<pxr::GfVec2f> uv_mesh, uv_result;
            pxr::VtArray<OpenMesh::Vec2f> flatxy;
            double energy;
            int fix_id1, fix_id2;
	};
}