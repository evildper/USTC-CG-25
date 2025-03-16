#include "freehand.h"
#include <iostream>
#include <imgui.h>

namespace USTC_CG
{
// Draw the line using ImGui
// Freehand::Freehand(float x1,float y1)
// {
//     x_list_.push_back(x1);
//     y_list_.push_back(y1);
// }
ImVec2 bezier(float t, float x0, float y0, float x1, float y1, float x2, float y2)
{
    float u = 1 - t;
    float xx = u * u * x0 + 2 * u * t * x1 + t * t * x2;
    float yy = u * u * y0 + 2 * u * t * y1 + t * t * y2;
    return ImVec2(xx, yy);
}

void Freehand::draw(const Config& config) const 
{
    //std::cout << "Total Points: " << x_list_.size() << std::endl;

    for (size_t i = 0; i + 2 < x_list_.size(); i=i+2)
    {
        for (float t = 0; t <= 1.0f; t += 0.01f)
        {
            ImVec2 p = bezier(t, x_list_[i], y_list_[i], x_list_[i + 1], y_list_[i + 1], x_list_[i + 2], y_list_[i + 2]);
            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->AddCircleFilled(ImVec2(p.x+config.bias[0],p.y+config.bias[1]), config.line_thickness/2, IM_COL32(config.line_color[0], config.line_color[1], config.line_color[2], config.line_color[3]));
        }
    }
}

void Freehand::add_control_point(float x, float y)
{
    x_list_.push_back(x);
    y_list_.push_back(y);
}
void Freehand::update(float x, float y)
{
    x_list_.push_back(x);
    y_list_.push_back(y);
}
}  // namespace USTC_CG