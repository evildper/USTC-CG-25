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
    // You can also compute pressure in this function 
    bool sifjasi=true;
    int n_particles = ps_.particles().size();
#pragma omp parallel for
    for (int i = 0; i < n_particles; i++) {
        auto p = ps_.particles()[i];
        //  ... necessary initialization of particle p's density here
        auto p_density = ps_.mass() * W_zero(ps_.h());
        auto density0 = ps_.density0();

        // Then traverse all neighbor fluid particles of p
        for (auto& q : p->neighbors()) {
            // ... compute the density contribution from q to p
            p_density += (ps_.mass()) * W(q->x() - p->x(), ps_.h());
        }
        p_density = std::max(p_density, density0);
        p->density() = p_density;

        p->pressure() = std::max(0.0, stiffness_ * (std::pow(p_density / density0, exponent_) - 1));
    }

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

    compute_density();
    compute_non_pressure_acceleration();
    compute_pressure_gradient_acceleration();        
    advect();


    TOC(step)
}
}  // namespace USTC_CG::node_sph_fluid