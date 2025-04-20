#include "path.h"
#include <algorithm>

#include <random>

#include "../surfaceInteraction.h"
USTC_CG_NAMESPACE_OPEN_SCOPE
using namespace pxr;

VtValue PathIntegrator::Li(const GfRay& ray, std::default_random_engine& random)
{
    std::uniform_real_distribution<float> uniform_dist(
        0.0f, 1.0f - std::numeric_limits<float>::epsilon());
    std::function<float()> uniform_float = std::bind(uniform_dist, random);

    auto color = EstimateOutGoingRadiance(ray, uniform_float, 0);

    return VtValue(GfVec3f(color[0], color[1], color[2]));
}
/*
 * TODO: You need to complete this function to achieve the estimate of the outgoing Radiance
 * */
pxr::GfVec3f ClampVec3f(const pxr::GfVec3f& v, float minVal, float maxVal) {
    return pxr::GfVec3f(
        std::max(minVal, std::min(maxVal, v[0])),
        std::max(minVal, std::min(maxVal, v[1])),
        std::max(minVal, std::min(maxVal, v[2]))
    );
}
GfVec3f PathIntegrator::EstimateOutGoingRadiance(
    const GfRay& ray,
    const std::function<float()>& uniform_float,
    int recursion_depth)
{
    if (recursion_depth >= 50) {
        return {};
    }

    SurfaceInteraction si;
    if (!Intersect(ray, si)) {
        if (recursion_depth == 0) {
            return IntersectDomeLight(ray);
        }

        return GfVec3f{ 0, 0, 0 };
    }

    // This can be customized : Do we want to see the lights? (Other than dome
    // lights?)

    


    // Flip the normal if opposite
    if (GfDot(si.shadingNormal, ray.GetDirection()) > 0) {
        si.flipNormal();
        si.PrepareTransforms();
    }

    GfVec3f color{ 0 };
    GfVec3f directLight = EstimateDirectLight(si, uniform_float);
    const bool see_lights = true;
    if (see_lights && recursion_depth == 0) {
        GfVec3f pos = {NAN, NAN, NAN};
        auto light = IntersectLights(ray, pos);
        if (!isnan(pos[0]) && pos.GetLength() < 100.0)
            return light;
    }
    // HW7_TODO: Estimate global lighting here.
    float P_RR = 0.99f;
    if (uniform_float() > P_RR)
        return { 0.0f, 0.0f, 0.0f };

    GfVec3f globalLight = GfVec3f{0.f};
    
    GfVec3f dir;
    float pdf;
    si.Sample(dir, pdf, uniform_float);
    
    float eps = 1e-5;
    
    if(pdf < eps){
        return { 0.0f, 0.0f, 0.0f };
    }
    
    GfRay nextRay = GfRay(si.position + 0.0001f * si.geometricNormal, dir);

    GfVec3f wo = GfVec3f(-ray.GetDirection()[0], -ray.GetDirection()[1], -ray.GetDirection()[2]);
    GfVec3f wi = dir.GetNormalized();
    wo = wo.GetNormalized();

    float cosVal = GfDot(si.shadingNormal.GetNormalized(), dir.GetNormalized());

    GfVec3f fr = si.Eval(wi);

    globalLight =
        GfCompMult(EstimateOutGoingRadiance(nextRay, uniform_float, recursion_depth + 1), fr) *
        cosVal / pdf / P_RR;
    color = directLight + globalLight;
    //color = ClampVec3f(color, 0.0f, 1.0f);

    return color;
}

USTC_CG_NAMESPACE_CLOSE_SCOPE
