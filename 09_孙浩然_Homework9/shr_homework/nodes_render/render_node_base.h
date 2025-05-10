#pragma once
#include "../camera.h"
#include "../geometries/mesh.h"
#include "../global_payload_gl.hpp"
#include "../material.h"
#include "GL/GLResources.hpp"
#include "nodes/core/def/node_def.hpp"

USTC_CG_NAMESPACE_OPEN_SCOPE
inline void render_node_type_base(NodeTypeInfo* ntype)
{
    ntype->color[0] = 114 / 255.f;
    ntype->color[1] = 94 / 255.f;
    ntype->color[2] = 29 / 255.f;
    ntype->color[3] = 1.0f;
}
#define global_payload      params.get_global_payload<RenderGlobalPayloadGL&>()
#define instance_collection global_payload.InstanceCollection
inline ResourceAllocator& get_resource_allocator(ExeParams& params)
{
    return global_payload.resource_allocator;
}

#define resource_allocator get_resource_allocator(params)
#define shader_factory     get_shader_factory(params)
inline Hd_USTC_CG_Camera* get_free_camera(ExeParams& params)
{
    auto& cameras = global_payload.get_cameras();

    Hd_USTC_CG_Camera* free_camera;
    for (auto camera : cameras) {
        if (camera->GetId() != SdfPath::EmptyPath()) {
            free_camera = camera;
            break;
        }
    }
    return free_camera;
}

#define materials global_payload.get_materials()
#define meshes    global_payload.get_meshes()
#define lights    global_payload.get_lights()

USTC_CG_NAMESPACE_CLOSE_SCOPE
