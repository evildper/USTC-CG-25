//
// Copyright 2017 Pixar
//
// Licensed under the Apache License, Version 2.0 (the "Apache License")
// with the following modification; you may not use this file except in
// compliance with the Apache License and the following modification to it:
// Section 6. Trademarks. is deleted and replaced with:
//
// 6. Trademarks. This License does not grant permission to use the trade
//    names, trademarks, service marks, or product names of the Licensor
//    and its affiliates, except as required to comply with Section 4(c) of
//    the License and to reproduce the content of the NOTICE file.
//
// You may obtain a copy of the Apache License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the Apache License with the above modification is
// distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied. See the Apache License for the specific
// language governing permissions and limitations under the Apache License.
//
#ifndef PXR_IMAGING_PLUGIN_HD_EMBREE_MESH_SAMPLERS_H
#define PXR_IMAGING_PLUGIN_HD_EMBREE_MESH_SAMPLERS_H
#include <embree4/rtcore.h>

#include <bitset>

#include "../api.h"
#include "../sampler.h"
#include "pxr/base/vt/types.h"
#include "pxr/imaging/hd/meshUtil.h"
#include "pxr/pxr.h"

USTC_CG_NAMESPACE_OPEN_SCOPE
using namespace pxr;
/// \class Hd_USTC_CG_RTCBufferAllocator
///
/// Utility class to track which embree user vertex buffers are currently
/// in use.
class Hd_USTC_CG_RTCBufferAllocator {
   public:
    /// Constructor. By default, set everything to unallocated.
    Hd_USTC_CG_RTCBufferAllocator() : _bitset(0)
    {
    }

    /// Allocate a buffer by finding the first clear bit, using that as
    /// the buffer number, and setting the bit to mark it as used.
    /// \return An unused RTC user vertex buffer id, or -1 on failure.
    int Allocate();

    /// Free a buffer by clearing its bit.
    /// \param buffer The buffer to mark as unused.
    void Free(int buffer);

    /// Query how many buffers are currently in user for this geometry
    unsigned int NumBuffers();

    /// As of Embree3 the number of buffers was greatly increased
    /// however the maximum is only defined locally to the library
    /// as of v3.4.0 this was the number.
    static constexpr int PXR_MAX_USER_VERTEX_BUFFERS = 16;

   private:
    std::bitset<PXR_MAX_USER_VERTEX_BUFFERS> _bitset;
};

// ----------------------------------------------------------------------
// The classes below implement the Hd_USTC_CG_PrimvarSampler interface for
// the different interpolation modes that hydra supports. In some cases,
// implementations are broken out by geometry type (e.g. triangles vs
// subdiv).

/// \class Hd_USTC_CG_ConstantSampler
///
/// This class implements the Hd_USTC_CG_PrimvarSampler interface for primvars
/// with "constant" interpolation mode. This means that the buffer only has
/// one item, which should be returned for any (element, u, v) tuple.
class Hd_USTC_CG_ConstantSampler : public Hd_USTC_CG_PrimvarSampler {
   public:
    /// Constructor.
    /// \param name The name of the primvar.
    /// \param value The buffer data for the primvar.
    Hd_USTC_CG_ConstantSampler(const TfToken& name, const VtValue& value)
        : _buffer(name, value),
          _sampler(_buffer)
    {
    }

    /// Sample the primvar at an (element, u, v) location. For constant
    /// primvars, the buffer only contains one item, so we always return
    /// that item.
    /// \param element The element index to sample.
    /// \param u The u coordinate to sample.
    /// \param v The v coordinate to sample.
    /// \param value The memory to write the value to (only written on success).
    /// \param dataType The HdTupleType describing element values.
    /// \return True if the value was successfully sampled.
    bool Sample(
        unsigned int element,
        float u,
        float v,
        void* value,
        HdTupleType dataType) const override;

   private:
    const HdVtBufferSource _buffer;
    const Hd_USTC_CG_BufferSampler _sampler;
};

/// \class Hd_USTC_CG_UniformSampler
///
/// This class implements the Hd_USTC_CG_PrimvarSampler interface for primvars
/// with "uniform" interpolation mode. This means that the buffer has one
/// item per authored face. For unrefined meshes, Hd_USTC_CG_ will convert
/// mesh polygons to triangles, so this class optionally takes an array
/// called "primitiveParams" which maps from the face index embree reports
/// to the original authored face in the scene data. If primitiveParams is not
/// provided, this translation step is skipped.
class Hd_USTC_CG_UniformSampler : public Hd_USTC_CG_PrimvarSampler {
   public:
    /// Constructor.
    /// \param name The name of the primvar.
    /// \param value The buffer data for the primvar.
    /// \param primitiveParams A mapping from geometry face index to authored
    ///                        face index.
    Hd_USTC_CG_UniformSampler(
        const TfToken& name,
        const VtValue& value,
        const VtIntArray& primitiveParams)
        : _buffer(name, value),
          _sampler(_buffer),
          _primitiveParams(primitiveParams)
    {
    }

    /// Constructor.
    /// \param name The name of the primvar.
    /// \param value The buffer data for the primvar.
    Hd_USTC_CG_UniformSampler(const TfToken& name, const VtValue& value)
        : _buffer(name, value),
          _sampler(_buffer)
    {
    }

    /// Sample the primvar at an (element, u, v) location. For uniform
    /// primvars, optionally look up the authored face index in
    /// _primitiveParams[element] (which is stored encoded); then return
    /// _buffer[element].
    ///
    /// \param element The element index to sample.
    /// \param u The u coordinate to sample.
    /// \param v The v coordinate to sample.
    /// \param value The memory to write the value to (only written on success).
    /// \param dataType The HdTupleType describing element values.
    /// \return True if the value was successfully sampled.
    bool Sample(
        unsigned int element,
        float u,
        float v,
        void* value,
        HdTupleType dataType) const override;

   private:
    const HdVtBufferSource _buffer;
    const Hd_USTC_CG_BufferSampler _sampler;
    const VtIntArray _primitiveParams;
};

/// \class Hd_USTC_CG_TriangleVertexSampler
///
/// This class implements the Hd_USTC_CG_PrimvarSampler interface for primvars
/// on triangle meshes with "vertex" or "varying" interpolation modes. This
/// means the buffer has one item per vertex, and the result of sampling is a
/// barycentric interpolation of the hit face vertices. This class
/// requires the triangulated mesh topology, to map from the triangle index
/// (in "element") to the triangle vertices.
class Hd_USTC_CG_TriangleVertexSampler : public Hd_USTC_CG_PrimvarSampler {
   public:
    /// Constructor.
    /// \param name The name of the primvar.
    /// \param value The buffer data for the primvar.
    /// \param indices A map from triangle index to vertex indices in the
    ///                triangulated geometry.
    Hd_USTC_CG_TriangleVertexSampler(
        const TfToken& name,
        const VtValue& value,
        const VtVec3iArray& indices)
        : _buffer(name, value),
          _sampler(_buffer),
          _indices(indices)
    {
    }

    /// Sample the primvar at an (element, u, v) location. For vertex primvars,
    /// the vertex indices of the triangle are stored in _indices[element][0-2].
    /// After fetching the primvar value for each of the three vertices,
    /// they are interpolated as follows, per Embree specification:
    /// t_uv = (1-u-v)*t0 + u*t1 + v*t2
    ///
    /// \param element The element index to sample.
    /// \param u The u coordinate to sample.
    /// \param v The v coordinate to sample.
    /// \param value The memory to write the value to (only written on success).
    /// \param dataType The HdTupleType describing element values.
    /// \return True if the value was successfully sampled.
    bool Sample(
        unsigned int element,
        float u,
        float v,
        void* value,
        HdTupleType dataType) const override;

   private:
    const HdVtBufferSource _buffer;
    const Hd_USTC_CG_BufferSampler _sampler;
    const VtVec3iArray _indices;
};

/// \class Hd_USTC_CG_TriangleFaceVaryingSampler
///
/// This class implements the Hd_USTC_CG_PrimvarSampler interface for primvars
/// on triangle meshes with "face-varying" interpolation modes. This means that
/// each vertex of each face gets its own buffer item: vertex 0 as part of
/// face 0 might have value 1.0f, but vertex 0 as part of face 1 might have
/// value 2.0f. The primvar's memory layout is grouped by face, with one item
/// per vertex.
///
/// Concretely, a cube with 8 vertices would have 24 items
/// (6 faces * 4 vertices) in a face-varying primvar, and the index of the
/// item for face 2, vertex 3, would be (2 * 4 + 3) = 11.
///
/// Face-varying primvars are provided to the sampler un-triangulated, but
/// the size of the buffer is tied to the size of the topology, so
/// this class triangulates the input buffer before sampling.
class Hd_USTC_CG_TriangleFaceVaryingSampler : public Hd_USTC_CG_PrimvarSampler {
   public:
    /// Constructor. Triangulates the provided buffer data.
    /// \param name The name of the primvar.
    /// \param value The buffer data for the primvar.
    /// \param meshUtil An HdMeshUtil instance that knows how to triangulate
    ///                 the input buffer data.
    Hd_USTC_CG_TriangleFaceVaryingSampler(
        const TfToken& name,
        const VtValue& value,
        HdMeshUtil& meshUtil)
        : _buffer(name, _Triangulate(name, value, meshUtil)),
          _sampler(_buffer)
    {
    }

    /// Sample the primvar at an (element, u, v) location. For face varying
    /// primvars, the vertex indices are simply (element * 3 + 0->2), since
    /// all faces are triangles. After fetching the primvar value for each of
    /// the three vertices, they are interpolated as follows, per Embree
    /// specification:
    /// t_uv = (1-u-v)*t0 + u*t1 + v*t2
    ///
    /// \param element The element index to sample.
    /// \param u The u coordinate to sample.
    /// \param v The v coordinate to sample.
    /// \param value The memory to write the value to (only written on success).
    /// \param dataType The HdTupleType describing element values.
    /// \return True if the value was successfully sampled.
    bool Sample(
        unsigned int element,
        float u,
        float v,
        void* value,
        HdTupleType dataType) const override;

   private:
    const HdVtBufferSource _buffer;
    const Hd_USTC_CG_BufferSampler _sampler;

    // Pass the "value" parameter through HdMeshUtils'
    // ComputeTriangulatedFaceVaryingPrimvar(), which adjusts the primvar
    // buffer data for the triangulated topology. HdMeshUtil is provided
    // the source topology at construction time, so this class doesn't need
    // to provide it.
    static VtValue _Triangulate(
        const TfToken& name,
        const VtValue& value,
        HdMeshUtil& meshUtil);
};

/// \class Hd_USTC_CG_SubdivVertexSampler
///
/// This class implements the Hd_USTC_CG_PrimvarSampler interface for primvars
/// on subdiv meshes with "vertex" interpolation mode. This means the buffer has
/// one item per vertex, and the result of sampling is a reconstruction using
/// the subdivision scheme basis weights. It uses embree's user vertex buffers
/// and rtcInterpolate API to accomplish the sampling.
class Hd_USTC_CG_SubdivVertexSampler : public Hd_USTC_CG_PrimvarSampler {
   public:
    /// Constructor. Allocates an embree user vertex buffer, and uploads
    /// the primvar data. Only float-based types (float, GfVec3f, GfMatrix4f)
    /// are allowed, and embree has an exhaustible number of user vertex
    /// buffers (16 at last count).
    ///
    /// \param name The name of the primvar.
    /// \param value The buffer data for the primvar.
    /// \param meshScene The owning mesh's embree prototype scene.
    /// \param meshId The owning mesh's geometry id in the prototype scene.
    /// \param allocator A mesh-global object that tracks buffer usage.
    Hd_USTC_CG_SubdivVertexSampler(
        const TfToken& name,
        const VtValue& value,
        RTCScene meshScene,
        unsigned meshId,
        Hd_USTC_CG_RTCBufferAllocator* allocator);

    /// Destructor. Frees the embree user vertex buffer.
    ~Hd_USTC_CG_SubdivVertexSampler() override;

    /// Sample the primvar at an (element, u, v) location. This implementation
    /// delegates to rtcInterpolate(). Only float-based types (float, GfVec3f,
    /// GfMatrix4f) are allowed.
    ///
    /// \param element The element index to sample.
    /// \param u The u coordinate to sample.
    /// \param v The v coordinate to sample.
    /// \param value The memory to write the value to (only written on success).
    /// \param dataType The HdTupleType describing element values.
    /// \return True if the value was successfully sampled.
    bool Sample(
        unsigned int element,
        float u,
        float v,
        void* value,
        HdTupleType dataType) const override;

   private:
    int _embreeBufferId;
    const HdVtBufferSource _buffer;
    RTCScene _meshScene;
    unsigned _meshId;
    Hd_USTC_CG_RTCBufferAllocator* _allocator;
};

USTC_CG_NAMESPACE_CLOSE_SCOPE

#endif  // PXR_IMAGING_PLUGIN_HD_EMBREE_MESH_SAMPLERS_H
