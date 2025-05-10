#include "wcsph.h"
#include <iostream>
using namespace Eigen;

namespace USTC_CG::sph_fluid {

WCSPH::WCSPH(const MatrixXd& X, const Vector3d& box_min, const Vector3d& box_max)
    : SPHBase(X, box_min, box_max)
{
}

void WCSPH::compute_density()
{
	// -------------------------------------------------------------
	// (HW TODO) Implement the density computation
    this->SPHBase::compute_density();
    for (auto& p : ps().particles()) {
        p->pressure_ = std::max(0.0, stiffness_ * (pow(p->density() / ps_.density0(), exponent_) - 1));
	}
    // You can also compute pressure in this function 
	// -------------------------------------------------------------
}

void WCSPH::step()
{
    TIC(step)
    // -------------------------------------------------------------
    // (HW TODO) Follow the instruction in documents and PPT, 
    // implement the pipeline of fluid simulation 
    // -------------------------------------------------------------

	// Search neighbors, compute density, advect, solve pressure acceleration, etc. 
    ps_.assign_particles_to_cells();
    ps_.search_neighbors();

    this->compute_density();
    this->compute_non_pressure_acceleration();
    this->compute_pressure_gradient_acceleration();

    this->advect();


    TOC(step)
}
}  // namespace USTC_CG::node_sph_fluid