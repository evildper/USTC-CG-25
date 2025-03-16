#include <istream>
#include "Ellipse.h"
#include <cmath>
#include <imgui.h>
float ImVec2Distance(const ImVec2& a, const ImVec2& b)
{
    float dx = b.x - a.x;
    float dy = b.y - a.y;
    return std::sqrt(dx * dx + dy * dy);
}
namespace USTC_CG
{
// Draw the line using ImGui
void Ellipse::draw(const Config& config) const
{   
    

    ImDrawList* draw_list = ImGui::GetWindowDrawList();
    ImVec2 a=ImVec2(
        config.bias[0] + start_point_x_, config.bias[1] + start_point_y_);
    ImVec2 b=ImVec2(config.bias[0] + end_point_x_, config.bias[1] + end_point_y_);

    float rot = atan2(b.y - a.y, b.x - a.x);
    float distance=ImVec2Distance(a,b)/2;

    
    double c=0.6;
    draw_list->AddEllipse(
        ImVec2((a.x+b.x)/2,(a.y+b.y)/2),
        ImVec2(distance,0.6f*distance),
        IM_COL32(
            config.line_color[0],
            config.line_color[1],
            config.line_color[2],
            config.line_color[3]),
            rot,
            ImDrawFlags_None,
        config.line_thickness);
}

void Ellipse::update(float x, float y)
{
    end_point_x_ = x;
    end_point_y_ = y;
}
}  // namespace USTC_CG