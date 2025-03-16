#include "polygon.h"
#include "imgui.h"
#include <algorithm>
#include <vector>
#include <algorithm>
#include <iostream>
namespace USTC_CG
{


// 更新当前多边形（鼠标移动时调用）
void Polygon::add_control_point(float x, float y)
{
    x_list_.push_back(x);
    y_list_.push_back(y);
}

// 绘制多边形
void Polygon::draw(const Config& config) const
{
    if (x_list_.size() < 2 && y_list_.size() < 2)
        return; // 至少两个点才能画出一条边
        

    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    
    for (size_t i = 0; i < x_list_.size() - 1; i++)
    {
        draw_list->AddLine(
            ImVec2(config.bias[0] + x_list_[i], config.bias[1] + y_list_[i]),
            ImVec2(config.bias[0] + x_list_[i + 1], config.bias[1] + y_list_[i + 1]),
            IM_COL32(config.line_color[0], config.line_color[1], config.line_color[2], config.line_color[3]),
            config.line_thickness);
    }

    // 如果多边形闭合，则连接首尾
    if (x_list_.size() > 2)
    {
        draw_list->AddLine(
            ImVec2(config.bias[0] + x_list_.back(), config.bias[1] + y_list_.back()),
            ImVec2(config.bias[0] + x_list_.front(), config.bias[1] + y_list_.front()),
            IM_COL32(config.line_color[0], config.line_color[1], config.line_color[2], config.line_color[3]),
            config.line_thickness);
    }
}
void Polygon::update(float x,float y){
    return;
}
std::vector<std::pair<int, int>> Polygon::get_interior_pixels() const
{
    std::vector<std::pair<int, int>> int_pixels;
    if (x_list_.size() < 3) 
        return int_pixels;

    // 计算多边形的边界框
    float min_x = *std::min_element(x_list_.begin(), x_list_.end());
    float max_x = *std::max_element(x_list_.begin(), x_list_.end());
    float min_y = *std::min_element(y_list_.begin(), y_list_.end());
    float max_y = *std::max_element(y_list_.begin(), y_list_.end());

    // 转换为整数像素范围
    int start_x = static_cast<int>(std::floor(min_x));
    int end_x = static_cast<int>(std::ceil(max_x));
    int start_y = static_cast<int>(std::floor(min_y));
    int end_y = static_cast<int>(std::ceil(max_y));

    // 遍历边界框范围内的每个像素点
    for (int y = start_y; y <= end_y; ++y)
    {
        std::vector<float> intersections;

        // 计算该扫描线与多边形边的交点
        for (size_t i = 0; i < x_list_.size(); ++i)
        {
            size_t j = (i + 1) % x_list_.size();  // 邻接点
            float x1 = x_list_[i], y1 = y_list_[i];
            float x2 = x_list_[j], y2 = y_list_[j];

            // 计算扫描线与 (x1, y1) - (x2, y2) 线段的交点
            if ((y1 <= y && y2 > y) || (y2 <= y && y1 > y)) 
            {
                float x_intersection = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
                intersections.push_back(x_intersection);
            }
        }

        // 对交点排序
        std::sort(intersections.begin(), intersections.end());

        // 按奇偶规则填充
        for (size_t k = 0; k + 1 < intersections.size(); k += 2) 
        {
            int x_start = static_cast<int>(std::ceil(intersections[k]));
            int x_end = static_cast<int>(std::floor(intersections[k + 1]));

            for (int x = x_start; x <= x_end; ++x) 
            {
                int_pixels.emplace_back(x, y);
            }
        }
        std::cout<<"the size of the mask is "<<int_pixels.size()<<std::endl;
    }
    return int_pixels;
}

}  // namespace USTC_CG