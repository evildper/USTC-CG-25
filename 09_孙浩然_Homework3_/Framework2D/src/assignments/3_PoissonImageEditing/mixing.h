#pragma once
#include "source_image_widget.h"
#include "common/image_widget.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
namespace USTC_CG
{
    class MixingClone
    {
       public:
        bool judge_A_precomputed(); // 判断 A 是否已经分解
        
        std::shared_ptr<Image> solve(); // 给外部调用的接口，求解 Poisson 方程组，返回一个 Seamless Clone 的结果图像（和背景图像一样大，替换了选中区域）
        void set_size(int width,int heigth); // 返回选中区域在背景图像中的位置
        void set_offset(int x,int y); // 设置选中区域的大小
        Eigen::SparseMatrix<double> A;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;  // 预分解求解器
        bool A_precomputed_ = false;  // 标记 A 是否已经分解
        void precompute_A();
        void compute_bounding_box();
        void set_right_to_mouse(int a,int b); // 预分解 A 矩阵
        void set_images(std::shared_ptr<Image> src,std::shared_ptr<Image> tar,std::shared_ptr<Image> mask); // 设置源图像、背景图像和选中区域
       private:
        // 注意使用指针，避免额外的复制操作
        
        
        std::shared_ptr<Image> src_img_; // 源图像
        std::shared_ptr<Image> tar_img_; // 背景图像
        std::shared_ptr<Image> src_selected_mask_; // 选中区域（矩形情形可以无视）
        int offset_x_, offset_y_;        // 矩形区域在背景图像中的位置（例如，左上角的坐标）
        int width_, height_; 
        int right_to_mouse_x,right_to_mouse_y;              // 选中区域的宽度和高度
    };
}