#pragma once

#include "shape.h"
#include <vector>
namespace USTC_CG
{
class Polygon : public Shape
{
   public:
    Polygon()=default;
    Polygon(float x, float y) { add_control_point(x, y); }
    virtual ~Polygon() = default;

    void draw(const Config& config) const override;
    void update(float x, float y);
    void add_control_point(float x, float y);

   private:
    std::vector<float> x_list_, y_list_;
};
}
