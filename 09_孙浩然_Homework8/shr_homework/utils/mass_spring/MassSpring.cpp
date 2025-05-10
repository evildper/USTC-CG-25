#include "MassSpring.h"
#include <iostream>

namespace USTC_CG::mass_spring {
MassSpring::MassSpring(const Eigen::MatrixXd& X, const EdgeSet& E)
{
    this->X = this->init_X = X;
    this->vel = Eigen::MatrixXd::Zero(X.rows(), X.cols());
    this->E = E;

    std::cout << "number of edges: " << E.size() << std::endl;
    std::cout << "init mass spring" << std::endl;

    // Compute the rest pose edge length
    for (const auto& e : E) {
        Eigen::Vector3d x0 = X.row(e.first);
        Eigen::Vector3d x1 = X.row(e.second);
        this->E_rest_length.push_back((x0 - x1).norm());
    }

    // Initialize the mask for Dirichlet boundary condition
    dirichlet_bc_mask.resize(X.rows(), false);

    // (HW_TODO) Fix two vertices, feel free to modify this 
    unsigned n_fix = sqrt(X.rows());  // Here we assume the cloth is square
    dirichlet_bc_mask[0] = true;
    dirichlet_bc_mask[n_fix - 1] = true;
    // for (int i = 0; i < n_fix; i++)
    // {
    //     dirichlet_bc_mask[i] = true;
    // }



}

void MassSpring::step()

{   

    Eigen::Vector3d acceleration_ext = gravity + wind_ext_acc;

    unsigned n_vertices = X.rows();

    // The reason to not use 1.0 as mass per vertex: the cloth gets heavier as we increase the resolution
    double mass_per_vertex =
        mass / n_vertices; 

    //----------------------------------------------------
    // (HW Optional) Bonus part: Sphere collision
    Eigen::MatrixXd acceleration_collision =
        getSphereCollisionForce(sphere_center.cast<double>(), sphere_radius);
    //----------------------------------------------------

    if (time_integrator == IMPLICIT_EULER) {
        // Implicit Euler
        TIC(step)

        // (HW TODO)
       
        auto H_elastic = computeHessianSparse(stiffness);  // size = [nx3, nx3]
        // compute Y
        MatrixXd y = X + h * vel;
        y.rowwise() += h * h * acceleration_ext.transpose();
        // collision
        if (enable_sphere_collision) {
            y += h * h * acceleration_collision;
        }
        // compute grad
        auto grad = computeGrad(stiffness);
        grad = flatten(grad);

        // Solve Newton's search direction with linear solver
        X = flatten(X);
        y = flatten(y);
        MatrixXd b = -(mass_per_vertex / h / h * (X - y) + grad);
        for (int i = 0; i < n_vertices; i++) {
            if (dirichlet_bc_mask[i]) {
                for (int k = 0; k < 3; k++) {
                    b(3 * i + k) = 0;
                }
            }
        }
        Eigen::SparseLU<Eigen::SparseMatrix<double>> cg;
        cg.compute(H_elastic);
        MatrixXd s = cg.solve(b);
        X = unflatten(X);
        s = unflatten(s);
        // update X and vel
        vel = s / h;
        if (enable_damping)
            vel *= damping;
        X += s;
        TOC(step)
    }
    
    else if (time_integrator == SEMI_IMPLICIT_EULER) {
        TIC(step)
        // Semi-implicit Euler
        Eigen::MatrixXd acceleration = -computeGrad(stiffness) / mass_per_vertex;
        acceleration.rowwise() += acceleration_ext.transpose();

        // -----------------------------------------------
        // (HW Optional)
        if (enable_sphere_collision) {
            acceleration += acceleration_collision;
        }
        // -----------------------------------------------

        // (HW TODO): Implement semi-implicit Euler time integration

        Eigen::MatrixXd accel = acceleration / mass;
        vel += h * accel;
        if (enable_damping) {
            vel *= damping;
        }
        for (int i = 0; i < X.rows(); ++i) {
            if (dirichlet_bc_mask[i]) {
                // 对于固定点，速度和加速度都应该是零
                vel.row(i).setZero();  // 将速度设为零
                acceleration.row(i).setZero();  // 将外力设为零
            }
        }
        X += h * vel;
        // Update X and vel 
        TOC(step)   
    }
    else {
        std::cerr << "Unknown time integrator!" << std::endl;
        return;
    }
}

// There are different types of mass spring energy:
// For this homework we will adopt Prof. Huamin Wang's energy definition introduced in GAMES103
// course Lecture 2 E = 0.5 * stiffness * sum_{i=1}^{n} (||x_i - x_j|| - l)^2 There exist other
// types of energy definition, e.g., Prof. Minchen Li's energy definition
// https://www.cs.cmu.edu/~15769-f23/lec/3_Mass_Spring_Systems.pdf
double MassSpring::computeEnergy(double stiffness)
{
    double sum = 0.;
    unsigned i = 0;
    for (const auto& e : E) {
        auto diff = X.row(e.first) - X.row(e.second);
        auto l = E_rest_length[i];
        sum += 0.5 * stiffness * std::pow((diff.norm() - l), 2);
        i++;
    }
    return sum;
}

Eigen::MatrixXd MassSpring::computeGrad(double stiffness)
{
    Eigen::MatrixXd g = Eigen::MatrixXd::Zero(X.rows(), X.cols());
    unsigned i = 0;
    for (const auto& e : E) {
        // --------------------------------------------------
        // (HW TODO): Implement the gradient computation
        int idx0 = e.first;
        int idx1 = e.second;
        auto diff = X.row(idx0) - X.row(idx1);
        auto l = E_rest_length[i];
        auto norm = diff.norm();
        auto grad = stiffness * (norm - l) * diff / norm; // gradient of the energy w.r.t. X[idx0] and X[idx1]
        g.row(idx0) += grad;
        g.row(idx1) -= grad; // gradient of the energy w.r.t. X[idx0] and X[idx1]
        //const auto 也是有序的
        
        // --------------------------------------------------
        i++;
    }
    return g;
}

Eigen::SparseMatrix<double> MassSpring::computeHessianSparse(double stiffness)
{
    unsigned n_vertices = X.rows();
    Eigen::SparseMatrix<double> H(n_vertices * 3, n_vertices * 3); 

    bool enable_fix_points = !dirichlet_bc_mask.empty();  // ✅ 自动判断有没有固定点

    unsigned i = 0; 
    auto k = stiffness; 
    const auto I = Eigen::MatrixXd::Identity(3, 3); 
    std::vector<Eigen::Triplet<double>> tripletList;

    for (const auto& e : E) 
    {
        // --------------------------------------------------
        // (HW TODO): Implement the sparse version Hessian computation
        // Remember to consider fixed points
        // You can also consider positive definiteness here
        Eigen::MatrixXd x = X.row(e.first) - X.row(e.second);
        double l = x.norm();
        Eigen::Matrix3d H_e =
            stiffness * x.transpose() * x / l / l +
            stiffness * (1 - E_rest_length[i]/l) * (I - x.transpose() * x / l / l);
        // Update the blocks of the Hessian matrix
        int start_row = 3 * e.first;
        int start_col = 3 * e.second;

        // Update the diagonal 
        if (!dirichlet_bc_mask[e.first]) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; k++) {
                    tripletList.push_back(Eigen::Triplet<double>(
                        start_row + j, start_row + k, H_e(j, k)));
                    tripletList.push_back(Eigen::Triplet<double>(
                        start_row + j, start_col + k, - H_e(j, k)));
                }
            }
        }
        
        if (!dirichlet_bc_mask[e.second]) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; k++) {
                    tripletList.push_back(
                        Eigen::Triplet<double>(start_col + j, start_row + k, -H_e(j, k)));
                    tripletList.push_back(
                        Eigen::Triplet<double>(start_col + j, start_col + k, H_e(j, k)));
                }
            }
        }
        
        // --------------------------------------------------
        i++;
    }
    for (int i = 0; i < n_vertices; i++)
    {
        if (dirichlet_bc_mask[i]) {
            for (int k = 0; k < 3; k++) {
                tripletList.push_back(Eigen::Triplet<double>(3 * i + k, 3 * i + k, 1));
            }
        }
        else
        {
            for (int k = 0; k < 3; k++) {
                tripletList.push_back(Eigen::Triplet<double>(3 * i + k, 3 * i + k, mass / n_vertices / h / h));
            }
        }
    }
    H.setFromTriplets(tripletList.begin(), tripletList.end());
    // check SPD
    if (!checkSPD(H))
    {
        H.diagonal().array() += 1e-8;
    }
    H.makeCompressed();
    return H;
}





bool MassSpring::checkSPD(const Eigen::SparseMatrix<double>& A)
{
    // Eigen::SimplicialLDLT<SparseMatrix_d> ldlt(A);
    // return ldlt.info() == Eigen::Success;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    auto eigen_values = es.eigenvalues();
    return eigen_values.minCoeff() >= 1e-10;
}

void MassSpring::reset()
{
    std::cout << "reset" << std::endl;
    this->X = this->init_X;
    this->vel.setZero();
}

// ----------------------------------------------------------------------------------
// (HW Optional) Bonus part
Eigen::MatrixXd MassSpring::getSphereCollisionForce(Eigen::Vector3d center, double radius)
{
    Eigen::MatrixXd force = Eigen::MatrixXd::Zero(X.rows(), X.cols());
    for (int i = 0; i < X.rows(); i++) {
       // (HW Optional) Implement penalty-based force here 
       Eigen::Vector3d x = { X(i, 0), X(i, 1), X(i, 2) };
       double a = (x - center).norm();
       if (collision_scale_factor * radius - a > 0) {
           force.row(i) +=
               collision_penalty_k * (x - center) / a * (collision_scale_factor * radius - a);
       }
    }
    return force;
}
// ----------------------------------------------------------------------------------
 
bool MassSpring::set_dirichlet_bc_mask(const std::vector<bool>& mask)
{
	if (mask.size() == X.rows())
	{
		dirichlet_bc_mask = mask;
		return true;
	}
	else
		return false;
}

bool MassSpring::update_dirichlet_bc_vertices(const MatrixXd &control_vertices)
{
   for (int i = 0; i < dirichlet_bc_control_pair.size(); i++)
   {
       int idx = dirichlet_bc_control_pair[i].first;
	   int control_idx = dirichlet_bc_control_pair[i].second;
	   X.row(idx) = control_vertices.row(control_idx);
   }

   return true; 
}

bool MassSpring::init_dirichlet_bc_vertices_control_pair(const MatrixXd &control_vertices,
    const std::vector<bool>& control_mask)
{
    
	if (control_mask.size() != control_vertices.rows())
			return false; 

   // TODO: optimize this part from O(n) to O(1)
   // First, get selected_control_vertices
   std::vector<VectorXd> selected_control_vertices; 
   std::vector<int> selected_control_idx; 
   for (int i = 0; i < control_mask.size(); i++)
   {
       if (control_mask[i])
       {
			selected_control_vertices.push_back(control_vertices.row(i));
            selected_control_idx.push_back(i);
		}
   }

   // Then update mass spring fixed vertices 
   for (int i = 0; i < dirichlet_bc_mask.size(); i++)
   {
       if (dirichlet_bc_mask[i])
       {
           // O(n^2) nearest point search, can be optimized
           // -----------------------------------------
           int nearest_idx = 0;
           double nearst_dist = 1e6; 
           VectorXd X_i = X.row(i);
           for (int j = 0; j < selected_control_vertices.size(); j++)
           {
               double dist = (X_i - selected_control_vertices[j]).norm();
               if (dist < nearst_dist)
               {
				   nearst_dist = dist;
				   nearest_idx = j;
			   }
           }
           //-----------------------------------------
           
		   X.row(i) = selected_control_vertices[nearest_idx];
           dirichlet_bc_control_pair.push_back(std::make_pair(i, selected_control_idx[nearest_idx]));
	   }
   }

   return true; 
}

}  // namespace USTC_CG::node_mass_spring

