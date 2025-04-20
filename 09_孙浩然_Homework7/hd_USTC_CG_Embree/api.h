
#pragma once

#define USTC_CG_NAMESPACE_OPEN_SCOPE  namespace USTC_CG {
#define USTC_CG_NAMESPACE_CLOSE_SCOPE }

#if defined(_MSC_VER)
#define HD_USTC_CG_EMBREE_EXPORT   __declspec(dllexport)
#define HD_USTC_CG_EMBREE_IMPORT   __declspec(dllimport)
#define HD_USTC_CG_EMBREE_NOINLINE __declspec(noinline)
#define HD_USTC_CG_EMBREE_INLINE   __forceinline
#else
#define HD_USTC_CG_EMBREE_EXPORT __attribute__((visibility("default")))
#define HD_USTC_CG_EMBREE_IMPORT
#define HD_USTC_CG_EMBREE_NOINLINE __attribute__((noinline))
#define HD_USTC_CG_EMBREE_INLINE   __attribute__((always_inline)) inline
#endif

#if BUILD_HD_USTC_CG_EMBREE_MODULE
#define HD_USTC_CG_EMBREE_API    HD_USTC_CG_EMBREE_EXPORT
#define HD_USTC_CG_EMBREE_EXTERN extern
#else
#define HD_USTC_CG_EMBREE_API HD_USTC_CG_EMBREE_IMPORT
#if defined(_MSC_VER)
#define HD_USTC_CG_EMBREE_EXTERN
#else
#define HD_USTC_CG_EMBREE_EXTERN extern
#endif
#endif
