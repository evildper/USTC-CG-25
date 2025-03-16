#include "polygon.h"
#include "imgui.h"
#include <vector>
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

}  // namespace USTC_CG
