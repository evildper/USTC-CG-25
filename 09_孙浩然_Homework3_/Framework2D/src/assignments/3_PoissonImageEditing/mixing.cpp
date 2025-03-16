
#include "mixing.h"
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace USTC_CG {

void MixingClone::set_images(std::shared_ptr<Image> src, std::shared_ptr<Image> tar, std::shared_ptr<Image> mask) {
    src_img_ = src;
    tar_img_ = tar;
    src_selected_mask_ = mask;
}
void MixingClone::set_right_to_mouse(int a,int b){
    right_to_mouse_x = a;
    right_to_mouse_y = b;
}
void MixingClone::compute_bounding_box() {
    int W = src_img_->width(), H = src_img_->height();
    int x_min = W, x_max = 0, y_min = H, y_max = 0;
    
    bool found = false;
    for (int i = 0; i < W; i++) {
        for (int j = 0; j < H; j++) {
            if (src_selected_mask_->get_pixel(i, j)[0] == 255) {
                x_min = std::min(x_min, i);
                x_max = std::max(x_max, i);
                y_min = std::min(y_min, j);
                y_max = std::max(y_max, j);
                found = true;
            }
        }
    }
    
    if (!found) {
        std::cerr << "Error: No valid mask found!" << std::endl;
        return;
    }
    
    width_ = x_max - x_min + 1;
    height_ = y_max - y_min + 1;
    offset_x_ = x_min;
    offset_y_ = y_min;
}

std::shared_ptr<Image> MixingClone::solve() {
    //compute_bounding_box();
    int W = width_, H = height_;
    auto is_boundary = [&](int x, int y) {
        if (src_selected_mask_->get_pixel(x, y)[0] == 0) return false; // 不是掩码区域
        return 
            (x > 0 && src_selected_mask_->get_pixel(x - 1, y)[0] == 0) ||
            (x < W - 1 && src_selected_mask_->get_pixel(x + 1, y)[0] == 0) ||
            (y > 0 && src_selected_mask_->get_pixel(x, y - 1)[0] == 0) ||
            (y < H - 1 && src_selected_mask_->get_pixel(x, y + 1)[0] == 0);
    };
    std::unordered_map<int, int> mask_index_map;
    int idx_count = 0;

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int global_x = x + offset_x_;//good thing
            int global_y = y + offset_y_;
            if (src_selected_mask_->get_pixel(global_x, global_y)[0] == 255) {
                mask_index_map[y * W + x] = idx_count++;
            }
        }
    }
    
    std::vector<Eigen::VectorXd> B_channels(3, Eigen::VectorXd( idx_count ));

    for (auto [flat_idx, local_idx] : mask_index_map) {
        int x = flat_idx % W, y = flat_idx / W;

        if(mask_index_map.count((y - 1) * W + x)&&mask_index_map.count((y + 1) * W + x)&&mask_index_map.count(y * W + (x - 1))&&mask_index_map.count(y * W + (x + 1))){
            int global_x = x + offset_x_,global_y = y+ offset_y_;
            for(int c = 0;c < 3;c++)
            {
                double src_x = src_img_->get_pixel(x + 1 + offset_x_, y + offset_y_)[c] - src_img_->get_pixel(x + offset_x_, y + offset_y_)[c];
                double tar_x = tar_img_->get_pixel(x + 1 + offset_x_+right_to_mouse_x, y + offset_y_+right_to_mouse_y)[c] - tar_img_->get_pixel(x + offset_x_ +right_to_mouse_x, y + offset_y_+right_to_mouse_y)[c];

                double src_y = src_img_->get_pixel(x + offset_x_, y + 1 + offset_y_)[c] - src_img_->get_pixel(x + offset_x_, y + offset_y_)[c];
                double tar_y = tar_img_->get_pixel(x + offset_x_ +right_to_mouse_x, y + 1 + offset_y_ +right_to_mouse_y)[c] - tar_img_->get_pixel(x + offset_x_+right_to_mouse_x, y + offset_y_+right_to_mouse_y)[c];                  
                    // 选择梯度较大的方向，但保持方向
                double g_x = (std::abs(src_x) > std::abs(tar_x)) ? src_x : tar_x;
                double g_y = (std::abs(src_y) > std::abs(tar_y)) ? src_y : tar_y;

                    // 计算左、右、上、下梯度（用于散度计算）
                double g_x_left = (std::abs(
                                        src_img_->get_pixel(x + offset_x_, y + offset_y_)[c] -
                                        src_img_->get_pixel(x - 1 + offset_x_, y + offset_y_)[c]) >
                                    std::abs(
                                        tar_img_->get_pixel(x + offset_x_+right_to_mouse_x, y + offset_y_ +right_to_mouse_y)[c] -
                                        tar_img_->get_pixel(x - 1 + offset_x_+right_to_mouse_x, y + offset_y_+right_to_mouse_y)[c]))
                                        ? (src_img_->get_pixel(x + offset_x_, y + offset_y_)[c] -
                                        src_img_->get_pixel(x - 1 + offset_x_, y + offset_y_)[c])
                                        : (tar_img_->get_pixel(x + offset_x_+right_to_mouse_x, y + offset_y_ +right_to_mouse_y)[c] -
                                        tar_img_->get_pixel(x - 1 + offset_x_+right_to_mouse_x, y + offset_y_+right_to_mouse_y)[c]);

                double g_y_up = (std::abs(
                                        src_img_->get_pixel(x + offset_x_, y + offset_y_)[c] -
                                        src_img_->get_pixel(x + offset_x_, y - 1 + offset_y_)[c]) >
                                    std::abs(
                                        tar_img_->get_pixel(x + offset_x_ + right_to_mouse_x, y + offset_y_ +right_to_mouse_y)[c] -
                                        tar_img_->get_pixel(x + offset_x_+right_to_mouse_x, y - 1 + offset_y_+right_to_mouse_y)[c]))
                                        ? (src_img_->get_pixel(x + offset_x_, y + offset_y_)[c] -
                                        src_img_->get_pixel(x + offset_x_, y - 1 + offset_y_)[c])
                                        : (tar_img_->get_pixel(x + offset_x_ + right_to_mouse_x, y + offset_y_ +right_to_mouse_y)[c] -
                                        tar_img_->get_pixel(x + offset_x_+right_to_mouse_x, y - 1 + offset_y_+right_to_mouse_y)[c]);

                    // 计算散度 div(v)
                double div_x = g_x - g_x_left;
                double div_y = g_y - g_y_up;

                    // 计算最终 B 值
                B_channels[c](local_idx) = -div_x - div_y;//div is the opposite direction.
            }
        }
        else{
            for(int c=0;c<3;c++)
            {
                B_channels[c](local_idx) = tar_img_->get_pixel(x+offset_x_+right_to_mouse_x,y+offset_y_+right_to_mouse_y)[c];
            }
        }

    }
    
    std::vector<Eigen::VectorXd> X_channels(3);
    for (int c = 0; c < 3; ++c) {
        X_channels[c] = solver.solve(B_channels[c]);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Solver failed for channel " << c << std::endl;
            return nullptr;
        }
    }

    for (auto [flat_idx, local_idx] : mask_index_map) {
        int x = flat_idx % W, y = flat_idx / W;
            
        unsigned char r = static_cast<unsigned char>(std::clamp(X_channels[0](local_idx), 0.0, 255.0));
        unsigned char g = static_cast<unsigned char>(std::clamp(X_channels[1](local_idx), 0.0, 255.0));
        unsigned char b = static_cast<unsigned char>(std::clamp(X_channels[2](local_idx), 0.0, 255.0));
            
        int tar_x = x + offset_x_ + right_to_mouse_x;
        int tar_y = y + offset_y_ + right_to_mouse_y;
        if (tar_x >= 0 && tar_x < tar_img_->width() && tar_y >= 0 && tar_y < tar_img_->height()) {
            tar_img_->set_pixel(tar_x, tar_y, {r, g, b});
        }
        
    }
    
    return tar_img_;
}
void MixingClone::precompute_A() {
    compute_bounding_box();
    int W = width_, H = height_;
    auto is_boundary = [&](int x, int y) {
        if (src_selected_mask_->get_pixel(x, y)[0] == 0) return false; // 不是掩码区域
        return 
            (x > 0 && src_selected_mask_->get_pixel(x - 1, y)[0] == 0) ||
            (x < W - 1 && src_selected_mask_->get_pixel(x + 1, y)[0] == 0) ||
            (y > 0 && src_selected_mask_->get_pixel(x, y - 1)[0] == 0) ||
            (y < H - 1 && src_selected_mask_->get_pixel(x, y + 1)[0] == 0);
    };
    //Eigen::SparseMatrix<double> A(W * H, W * H);
    std::vector<Eigen::Triplet<double>> triplet_list;
    std::unordered_map<int, int> mask_index_map;
    int idx_count = 0;

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int global_x = x + offset_x_;
            int global_y = y + offset_y_;
            if (src_selected_mask_->get_pixel(global_x, global_y)[0] == 255) {
                mask_index_map[y * W + x] = idx_count;
                idx_count = idx_count+1;
            }
        }
    }
    std::cout<<"THE IDX= "<<idx_count<<std::endl;
    Eigen::SparseMatrix<double> A(idx_count ,idx_count);
    if (mask_index_map.empty()) {
        std::cerr << "Error: No valid mask pixels found!" << std::endl;
        return;
    }

    for (auto [flat_idx, local_idx] : mask_index_map) {
        int x = flat_idx % W, y = flat_idx / W;
        if(mask_index_map.count((y - 1) * W + x)&&mask_index_map.count((y + 1) * W + x)&&mask_index_map.count(y * W + (x - 1))&&mask_index_map.count(y * W + (x + 1))){
            triplet_list.emplace_back(local_idx, local_idx, 4.0);
            triplet_list.emplace_back(local_idx, mask_index_map[(y - 1) * W + x], -1.0);
            triplet_list.emplace_back(local_idx, mask_index_map[(y + 1) * W + x], -1.0);
            triplet_list.emplace_back(local_idx, mask_index_map[y * W + (x - 1)], -1.0);
            triplet_list.emplace_back(local_idx, mask_index_map[y * W + (x + 1)], -1.0);
        }
        else{
            triplet_list.emplace_back(local_idx, local_idx, 1.0);
        }

    }

    
    A.setFromTriplets(triplet_list.begin(), triplet_list.end());

    // std::cout << "A.rows() = " << A.rows() << ", A.cols() = " << A.cols() << ", A.nonZeros() = " << A.nonZeros() << std::endl;

    // if (A.rows() == 0 || A.nonZeros() == 0) {
    //     std::cerr << "Error: A matrix is empty!" << std::endl;
    //     return;
    // }

    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Error: Solver precompute failed!" << std::endl;
        return;
    }

    std::cout << "A matrix initialized successfully!" << std::endl;
    A_precomputed_ = true;
}
bool MixingClone::judge_A_precomputed() {
    return A_precomputed_;
}

}  // namespace USTC_CG

// #include "mixing.h"
// #include <iostream>
// #include <Eigen/Sparse>
// #include <Eigen/Dense>

// namespace USTC_CG {

// void MixingClone::set_size(int width, int height) {
//     width_ = width;
//     height_ = height;
// }

// void MixingClone::set_offset(int x, int y) {
//     offset_x_ = x;
//     offset_y_ = y;
// }

// void MixingClone::set_images(std::shared_ptr<Image> src, std::shared_ptr<Image> tar, std::shared_ptr<Image> mask) {
//     src_img_ = src;
//     tar_img_ = tar;
//     src_selected_mask_ = mask;
// }
// std::shared_ptr<Image> MixingClone::solve() {
//     int W = width_, H = height_;
//     if (!src_img_ || !tar_img_ || !src_selected_mask_ || W <= 0 || H <= 0) {
//         std::cerr << "Invalid image or size" << std::endl;
//         return nullptr;
//     }
//     std::pair<int,int> p(0,0),q(0,0);
//     for(int i=0;i<width_;i++){
//         for(int j=0;j<height_;j++){
//             std::vector<unsigned char> pixelValues = src_selected_mask_->get_pixel(i, j);
//             if(pixelValues[0]==255){
//                     p.first = i;
//                     p.second = j;
//                     i = width_;
//                     j = height_;
//             }
//         }
//     }
//     for(int i=width_-1;i>=0;i--){
//         for(int j=height_-1;j>=0;j--){
//             std::vector<unsigned char> pixelValues = src_selected_mask_->get_pixel(i, j);
//             if(pixelValues[0]==255){
//                     q.first = i;
//                     q.second = j;
//                     i = -1;
//                     j = -1;
//             }
//         }
//     }
//     auto is_boundary = [&](int x, int y) {
//         return x == 0 || x == W - 1 || y == 0 || y == H - 1;
//     };
//     std::vector<Eigen::VectorXd> B_channels(3, Eigen::VectorXd(W * H));

//     for (int y = 0; y < H; y++) {
//         for (int x = 0; x < W; x++) {
//             int idx = y * W + x;
//             if (is_boundary(x, y)) {
//                 // 直接使用目标图像的颜色
//                 for (int c = 0; c < 3; ++c) {
//                     B_channels[c](idx) = tar_img_->get_pixel(x + offset_x_, y + offset_y_)[c];
//                 }
//             } 
//             else {
//                 // Mixing 计算混合梯度
//                 for (int c = 0; c < 3; ++c) {
//                     double src_x = src_img_->get_pixel(x + 1 + p.first, y + p.second)[c] - src_img_->get_pixel(x + p.first, y + p.second)[c];
//                     double tar_x = tar_img_->get_pixel(x + 1 + offset_x_, y + offset_y_)[c] - tar_img_->get_pixel(x + offset_x_, y + offset_y_)[c];

//                     double src_y = src_img_->get_pixel(x + p.first, y + 1 + p.second)[c] - src_img_->get_pixel(x + p.first, y + p.second)[c];
//                     double tar_y = tar_img_->get_pixel(x + offset_x_, y + 1 + offset_y_)[c] - tar_img_->get_pixel(x + offset_x_, y + offset_y_)[c];
//                     // double alpha = 0.001;
//                     // double g_x = alpha * src_x + (1 - alpha) * tar_x;
//                     // double g_y = alpha * src_y + (1 - alpha) * tar_y;
//                     // // 计算左、右、上、下梯度（用于散度计算）
//                     // double g_x_left = alpha*std::abs(
//                     //                     src_img_->get_pixel(x + p.first, y + p.second)[c] - 
//                     //                     src_img_->get_pixel(x - 1 + p.first, y + p.second)[c]) +(1-alpha)*
//                     //                 std::abs(
//                     //                     tar_img_->get_pixel(x + offset_x_, y + offset_y_)[c] - 
//                     //                     tar_img_->get_pixel(x - 1 + offset_x_, y + offset_y_)[c]);
//                     // double g_y_up = alpha*std::abs(
//                     //                     src_img_->get_pixel(x + p.first, y + p.second)[c] - 
//                     //                     src_img_->get_pixel(x + p.first, y - 1 + p.second)[c]) +(1-alpha)*
//                     //                 std::abs(
//                     //                     tar_img_->get_pixel(x + offset_x_, y + offset_y_)[c] - 
//                     //                     tar_img_->get_pixel(x + offset_x_, y - 1 + offset_y_)[c]);
//                     // 计算散度 div(v)


                                        
//                     // 选择梯度较大的方向，但保持方向
//                     double g_x = (std::abs(src_x) > std::abs(tar_x)) ? src_x : tar_x;
//                     double g_y = (std::abs(src_y) > std::abs(tar_y)) ? src_y : tar_y;

//                     // 计算左、右、上、下梯度（用于散度计算）
//                     double g_x_left = (std::abs(
//                                         src_img_->get_pixel(x + p.first, y + p.second)[c] -
//                                         src_img_->get_pixel(x - 1 + p.first, y + p.second)[c]) >
//                                     std::abs(
//                                         tar_img_->get_pixel(x + offset_x_, y + offset_y_)[c] -
//                                         tar_img_->get_pixel(x - 1 + offset_x_, y + offset_y_)[c]))
//                                         ? (src_img_->get_pixel(x + p.first, y + p.second)[c] - src_img_->get_pixel(x - 1 + p.first, y + p.second)[c])
//                                         : (tar_img_->get_pixel(x + offset_x_, y + offset_y_)[c] - tar_img_->get_pixel(x - 1 + offset_x_, y + offset_y_)[c]);

//                     double g_y_up = (std::abs(
//                                         src_img_->get_pixel(x + p.first, y + p.second)[c] -
//                                         src_img_->get_pixel(x + p.first, y - 1 + p.second)[c]) >
//                                     std::abs(
//                                         tar_img_->get_pixel(x + offset_x_, y + offset_y_)[c] -
//                                         tar_img_->get_pixel(x + offset_x_, y - 1 + offset_y_)[c]))
//                                         ? (src_img_->get_pixel(x + p.first, y + p.second)[c] - src_img_->get_pixel(x + p.first, y - 1 + p.second)[c])
//                                         : (tar_img_->get_pixel(x + offset_x_, y + offset_y_)[c] - tar_img_->get_pixel(x + offset_x_, y - 1 + offset_y_)[c]);

//                     // 计算散度 div(v)
//                     double div_x = g_x - g_x_left;
//                     double div_y = g_y - g_y_up;

//                     // 计算最终 B 值
//                     B_channels[c](idx) = -div_x - div_y;//div is the opposite direction.
//                 }
//             }
                
            
//         }
//     }

//     // 使用预分解的 A 快速求解
//     std::vector<Eigen::VectorXd> X_channels(3);
//     for (int c = 0; c < 3; ++c) {
//         X_channels[c] = solver.solve(B_channels[c]);
//         if (solver.info() != Eigen::Success) {
//             std::cerr << "Solver failed for channel " << c << std::endl;
//             return nullptr;
//         }
//     }
//     for (int y = 0; y < H; ++y) {
//         for (int x = 0; x < W; ++x) {
//             int idx = y * W + x;
//             unsigned char r = static_cast<unsigned char>(std::clamp(X_channels[0](idx), 0.0, 255.0));
//             unsigned char g = static_cast<unsigned char>(std::clamp(X_channels[1](idx), 0.0, 255.0));
//             unsigned char b = static_cast<unsigned char>(std::clamp(X_channels[2](idx), 0.0, 255.0));
//             // 计算在目标图像上的正确位置
//             int tar_x = x + offset_x_;
//             int tar_y = y + offset_y_;
//             // 仅在目标图像范围内修改
//             if (tar_x >= 0 && tar_x < tar_img_->width() && tar_y >= 0 && tar_y < tar_img_->height()) {
//                 tar_img_->set_pixel(tar_x, tar_y, {r, g, b});
//             }
//         }//not a image class so it must use the get data founction to get the image data.
//     }
//     return tar_img_;
    
// }
// void MixingClone::precompute_A() {
//     int W = width_, H = height_;
//     Eigen::SparseMatrix<double> A(W * H, W * H);
//     std::vector<Eigen::Triplet<double>> triplet_list;
    
//     auto is_boundary = [&](int x, int y) {
//         return x == 0 || x == W - 1 || y == 0 || y == H - 1;
//     };
    
//     for (int y = 0; y < H; y++) {
//         for (int x = 0; x < W; x++) {
//             int idx = y * W + x;
    
//             if (is_boundary(x, y)) {
//                 // 处理边界点
//                 triplet_list.emplace_back(idx, idx, 1.0);
//             } else {
//                 // 处理内部点 (Poisson 5-point Laplacian)
//                 triplet_list.emplace_back(idx, idx, 4.0);
//                 triplet_list.emplace_back(idx, (y - 1) * W + x, -1.0);
//                 triplet_list.emplace_back(idx, (y + 1) * W + x, -1.0);
//                 triplet_list.emplace_back(idx, y * W + (x - 1), -1.0);
//                 triplet_list.emplace_back(idx, y * W + (x + 1), -1.0);
//             }
//         }
//     }
    
//     A.setFromTriplets(triplet_list.begin(), triplet_list.end());
    
//     // 预分解 A
//     solver.compute(A);
//     if (solver.info() != Eigen::Success) {
//         std::cerr << "Precompute A failed!" << std::endl;
//         return;
//     }
//     std::cout << "A matrix nonzero elements: " << A.nonZeros() << std::endl;

//     A_precomputed_ = true;
// }

// bool MixingClone::judge_A_precomputed() {
//     return A_precomputed_;
// }

// }  // namespace USTC_CG