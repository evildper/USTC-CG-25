#pragma once

#include "shape.h"
#include <vector>
namespace USTC_CG
{
class Freehand : public Shape
{
   public:
    Freehand() = default;
    Freehand(float x, float y) { add_control_point(x, y); }
    
    // Constructor to initialize a line with start and end coordinates


    virtual ~Freehand() = default;

    // Overrides draw function to implement line-specific drawing logic
    void draw(const Config& config) const override;

    // Overrides Shape's update function to adjust the end point during
    // interaction
    void add_control_point(float x, float y) override;
    void update(float x, float y) override;

   private:
   std::vector<float> x_list_, y_list_;
};
}  // namespace USTC_CG