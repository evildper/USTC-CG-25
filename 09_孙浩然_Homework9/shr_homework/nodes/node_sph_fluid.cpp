#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <iostream>
#include <memory>

#include "GCore/Components/MeshOperand.h"
#include "GCore/Components/PointsComponent.h"
#include "GCore/util_openmesh_bind.h"
#include "geom_node_base.h"

#include "sph_fluid/wcsph.h"
#include "sph_fluid/iisph.h"
#include "sph_fluid/sph_base.h"

struct SPHFluidStorage {
  constexpr static bool has_storage = false;
  std::shared_ptr<USTC_CG::sph_fluid::SPHBase> sph_base;
};

// ------------------------- helper functions -------------------------------

pxr::VtArray<pxr::GfVec3f> eigen_to_usd_vertices(const Eigen::MatrixXd& V)
{
  pxr::VtArray<pxr::GfVec3f> vertices;
  for (int i = 0; i < V.rows(); i++) {
    vertices.push_back(pxr::GfVec3f(V(i, 0), V(i, 1), V(i, 2)));
  }
  return vertices;
}

Eigen::MatrixXd usd_vertices_to_eigen(const pxr::VtArray<pxr::GfVec3f>& V)
{
  Eigen::MatrixXd vertices(V.size(), 3);
  for (int i = 0; i < V.size(); i++) {
    vertices.row(i) << V[i][0], V[i][1], V[i][2];
  }
  return vertices;
}

// --------------------------------------------------------------------------

NODE_DEF_OPEN_SCOPE
NODE_DECLARATION_FUNCTION(sph_fluid) {
  b.add_input<Geometry>("Points");
  // Simulation parameters 
  b.add_input<pxr::GfVec3f>("sim box min");
  b.add_input<pxr::GfVec3f>("sim box max");

  // general parameters
  b.add_input<float>("dt").default_val(0.01).min(0.0).max(0.5);
  b.add_input<float>("viscosity").default_val(0.03).min(0.0).max(0.5);
  b.add_input<float>("gravity").default_val(-9.8).min(-20.0).max(20.0);

  // WCSPH parameters 
  b.add_input<float>("stiffness").default_val(500).min(100).max(10000);
  b.add_input<float>("exponent").default_val(7).min(1).max(10);

  // --------- (HW Optional) if you implement IISPH, please uncomment the following lines ------------

  //b.add_input<float>("omega").default_val(0.5).min(0.).max(1.);
  //b.add_input<int>("max iter").default_val(20).min(0).max(1000);

  // -----------------------------------------------------------------------------------------------------------

  // Useful switches (0 or 1). You can add more if you like.
  b.add_input<int>("enable time profiling").default_val(0).min(0).max(1);
  b.add_input<int>("enable debug output").default_val(0).min(0).max(1);

  // Optional switches
  b.add_input<int>("enable IISPH").default_val(0).min(0).max(1);

  // Output 
  b.add_output<Geometry>("Points");
  b.add_output<pxr::VtArray<pxr::GfVec3f>>("Point Colors");
}

NODE_EXECUTION_FUNCTION(sph_fluid) {
  using namespace Eigen;
  using namespace USTC_CG::sph_fluid;

  auto& global_payload = params.get_global_payload<GeomPayload&>();
  auto current_time = global_payload.current_time;

  auto& storage = params.get_storage<SPHFluidStorage&>();
  auto& sph_base = storage.sph_base;

  // ----------------------------- Load and check simulation box area ------------------------------------------------
  auto sim_box_min = params.get_input<pxr::GfVec3f>("sim box min");
  auto sim_box_max = params.get_input<pxr::GfVec3f>("sim box max");

  if (sim_box_max[0] <= sim_box_min[0] || sim_box_max[1] <= sim_box_min[1] || sim_box_max[2] <= sim_box_min[2]) {
    throw std::runtime_error("Invalid simulation box.");
  }
  // --------------------------- Load particles -------------------------------------------
  auto geometry = params.get_input<Geometry>("Points");
  auto points = geometry.get_component<PointsComponent>();
  if (!points || points->get_vertices().size() == 0) {
    throw std::runtime_error("Invalid point set.");
  }
  //---------------------------------------------------------------------------------------

  if (current_time == 0 || !sph_base) {  
    if (sph_base != nullptr)
      sph_base.reset();

    // Create particles positions
    MatrixXd particle_pos = usd_vertices_to_eigen(points->get_vertices());
    // Create simulation box (two end points in the space)
    Vector3d box_min{ sim_box_min[0], sim_box_min[1], sim_box_min[2] };
    Vector3d box_max{ sim_box_max[0], sim_box_max[1], sim_box_max[2] };
    // Create solver
    bool enable_IISPH = params.get_input<int>("enable IISPH") == 1 ? true : false;
    if (enable_IISPH) {
      sph_base = std::make_shared<IISPH>(particle_pos, box_min, box_max);
    }
    else {
      sph_base = std::make_shared<WCSPH>(particle_pos, box_min, box_max);
    }

    sph_base->dt() = params.get_input<float>("dt");
    sph_base->viscosity() = params.get_input<float>("viscosity");
    sph_base->gravity() = { 0, 0, params.get_input<float>("gravity") };

    // Useful switches
    sph_base->enable_time_profiling = params.get_input<int>("enable time profiling") == 1 ? true : false;
    sph_base->enable_debug_output = params.get_input<int>("enable debug output") == 1 ? true : false;

    if (enable_IISPH) {
      // --------- (HW Optional) if you implement IISPH please uncomment the following lines -----------		

          //std::dynamic_pointer_cast<IISPH>(sph_base)->max_iter() = params.get_input<int>("max iter");
          //std::dynamic_pointer_cast<IISPH>(sph_base)->omega() = params.get_input<float>("omega");

      // --------------------------------------------------------------------------------------------------------
    }
    else {
      std::dynamic_pointer_cast<WCSPH>(sph_base)->stiffness() = params.get_input<float>("stiffness");
      std::dynamic_pointer_cast<WCSPH>(sph_base)->exponent() = params.get_input<float>("exponent");
    }

  }
  else  // otherwise, step forward the simulation
  {
    sph_base->step();
  }

  // ------------------------- construct necessary output ---------------
  auto vertices = sph_base->getX();
  points->set_vertices(eigen_to_usd_vertices(vertices));

  auto color = eigen_to_usd_vertices(sph_base->get_vel_color_jet());

  params.set_output("Point Colors", std::move(color));
  params.set_output("Points", std::move(geometry));
  // ----------------------------------------------------------------------------------------
  return true;
}

NODE_DECLARATION_UI(sph_fluid);
NODE_DEF_CLOSE_SCOPE
