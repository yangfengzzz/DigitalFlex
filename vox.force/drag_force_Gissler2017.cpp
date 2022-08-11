//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/drag_force_Gissler2017.h"

#include "vox.force/boundary_model_Bender2019.h"
#include "vox.force/simulation.h"

namespace vox::flex {

// C.Gissler, S.Band, A.Peer, M.Ihmsen, M.Teschner,
// "Approximate Air-Fluid Interactions for SPH,"
// VRIPHYS 2017

DragForce_Gissler2017::DragForce_Gissler2017(FluidModel *model) : DragBase(model) {}

DragForce_Gissler2017::~DragForce_Gissler2017(void) = default;

#ifdef USE_AVX

void DragForce_Gissler2017::step() {
    Simulation *sim = Simulation::getCurrent();
    const double supportRadius = sim->getSupportRadius();
    const double radius = sim->getValue<double>(Simulation::PARTICLE_RADIUS);
    const double diam = static_cast<double>(2.0) * radius;
    static const double pi = static_cast<double>(M_PI);
    const double rho_l = m_model->getDensity0();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    FluidModel *model = m_model;
    const unsigned int numParticles = m_model->numActiveParticles();
    if (numParticles == 0) return;

    const double fluidParticleVolume = m_model->getVolume(0);

    // Air velocity.
    const Vector3D va(0, 0, 0);

    const double L = cbrt(static_cast<double>(0.75) / pi) * diam;

    const double inv_td = static_cast<double>(0.5) * C_d * mu_l / (rho_l * L * L);
    const double td = static_cast<double>(1.0) / inv_td;
    double omegaSquare = C_k * sigma / (rho_l * L * L * L) - inv_td * inv_td;
    omegaSquare = std::max(omegaSquare, static_cast<double>(0.0));
    const double omega = sqrt(omegaSquare);

    // Equation (6)
    double val = td * td * omegaSquare;
    val = sqrt(val + static_cast<double>(1.0)) + td * omega;
    val = std::max(val, -static_cast<double>(0.5) * pi);
    val = std::min(val, static_cast<double>(0.5) * pi);
    const double t_max = -static_cast<double>(2.0) * (atan(val) - pi) / omega;

    // Equation (7)
    const double c_def =
            static_cast<double>(1.0) -
            exp(-t_max / td) * (cos(omega * t_max) + static_cast<double>(1.0) / (omega * td) * sin(omega * t_max));

    // Weber number without velocity
    const double We_i_wo_v = rho_a * L / sigma;

    // Equation (8)
    const double y_coeff = (C_F * We_i_wo_v * c_def) / (C_k * C_b);

    const double n_full = 38;
    const double n_full_23 = n_full * 2.0 / 3.0;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)numParticles; i++) {
            const Vector3D &vi = m_model->getVelocity(i);
            Vector3D v_i_rel = va - vi;
            const double vi_rel_square = v_i_rel.squaredNorm();
            const double vi_rel_norm = sqrt(vi_rel_square);
            const double We_i = We_i_wo_v * vi_rel_square;

            Vector3D v_i_rel_n = v_i_rel;
            if (vi_rel_norm <= 1.0e-6) continue;
            // Else.
            v_i_rel_n = v_i_rel_n * (1.0 / vi_rel_norm);

            // Equation (8)
            const double y_i_max = std::min(vi_rel_square * y_coeff, static_cast<double>(1.0));

            const double Re_i =
                    static_cast<double>(2.0) * std::max((rho_a * vi_rel_norm * L) / mu_a, static_cast<double>(0.1));

            double C_Di_sphere;
            if (Re_i <= 1000.0)
                C_Di_sphere = static_cast<double>(24.0) / Re_i *
                              (static_cast<double>(1.0) +
                               static_cast<double>(1.0 / 6.0) * pow(Re_i, static_cast<double>(2.0 / 3.0)));
            else
                C_Di_sphere = static_cast<double>(0.424);

            // Equation (9)
            const double C_Di_Liu = C_Di_sphere * (static_cast<double>(1.0) + static_cast<double>(2.632) * y_i_max);

            unsigned int numNeighbors = 0;
            for (unsigned int pid = 0; pid < nFluids; pid++)
                numNeighbors += sim->numberOfNeighbors(fluidModelIndex, pid, i);

            if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
                for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++)
                    numNeighbors += sim->numberOfNeighbors(fluidModelIndex, pid, i);
            }
            // if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
            //{
            //	forall_density_maps(
            //		// ToDo
            //	);
            // }
            else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019) {
                forall_volume_maps(numNeighbors += (unsigned int)(Vj / fluidParticleVolume) + 1;);
            }

            // Equation (10)
            double C_Di;
            const double factor_n = std::min(n_full_23, (double)numNeighbors) / n_full_23;
            if (numNeighbors == 0)
                C_Di = C_Di_Liu;
            else
                C_Di = (static_cast<double>(1.0) - factor_n) * C_Di_Liu + factor_n;

            // Equation (12)
            const double h1 = (L + C_b * L * y_i_max);
            const double A_i_droplet = pi * h1 * h1;

            // Equation (13)
            const double A_i_unoccluded = (static_cast<double>(1.0) - factor_n) * A_i_droplet + factor_n * diam * diam;

            //////////////////////////////////////////////////////////////////////////
            // Fluid
            //////////////////////////////////////////////////////////////////////////
            Scalarf8 max_v_x_avx(0.0f);
            const Vector3D &xi = m_model->getPosition(i);

            const Vector3f8 xi_avx(xi);
            const Vector3f8 v_i_rel_n_avx(v_i_rel_n);

            forall_fluid_neighbors_in_same_phase_avx(Vector3f8 xixj = xi_avx - xj_avx; xixj.normalize();
                                                     const Scalarf8 x_v = v_i_rel_n_avx.dot(xixj);
                                                     max_v_x_avx = blend(max_v_x_avx > x_v, max_v_x_avx, x_v);)

                    float max_v_x_[8];
            max_v_x_avx.store(max_v_x_);
            double max_v_x = 0.0;
            for (int k = 0; k < 8; k++) max_v_x = std::max(max_v_x, max_v_x_[k]);

            // Equation (15)
            const double w_i = std::max(static_cast<double>(0.0),
                                        std::min(static_cast<double>(1.0), static_cast<double>(1.0) - max_v_x));

            // Equation (14)
            const double A_i = w_i * A_i_unoccluded;

            // Drag force. Additionally dividing by mass to get acceleration.
            Vector3D &ai = m_model->getAcceleration(i);
            ai += m_dragCoefficient * static_cast<double>(0.5) / m_model->getMass(i) * rho_a * (v_i_rel * vi_rel_norm) *
                  C_Di * A_i;
        }
    }
}

#else

void DragForce_Gissler2017::step() {
    Simulation *sim = Simulation::getCurrent();
    const auto radius = sim->getValue<double>(Simulation::PARTICLE_RADIUS);
    const double diam = static_cast<double>(2.0) * radius;
    static const auto pi = static_cast<double>(M_PI);
    const double rho_l = m_model->getDensity0();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    const unsigned int nFluids = sim->numberOfFluidModels();
    const unsigned int nBoundaries = sim->numberOfBoundaryModels();
    FluidModel *model = m_model;
    const unsigned int numParticles = m_model->numActiveParticles();
    if (numParticles == 0) return;

    const double fluidParticleVolume = m_model->getVolume(0);

    // Air velocity.
    const Vector3D va(0, 0, 0);

    const double L = cbrt(static_cast<double>(0.75) / pi) * diam;

    const double inv_td = static_cast<double>(0.5) * C_d * mu_l / (rho_l * L * L);
    const double td = static_cast<double>(1.0) / inv_td;
    double omegaSquare = C_k * sigma / (rho_l * L * L * L) - inv_td * inv_td;
    omegaSquare = std::max(omegaSquare, static_cast<double>(0.0));
    const double omega = sqrt(omegaSquare);

    // Equation (6)
    double val = td * td * omegaSquare;
    val = sqrt(val + static_cast<double>(1.0)) + td * omega;
    val = std::max(val, -static_cast<double>(0.5) * pi);
    val = std::min(val, static_cast<double>(0.5) * pi);
    const double t_max = -static_cast<double>(2.0) * (atan(val) - pi) / omega;

    // Equation (7)
    const double c_def =
            static_cast<double>(1.0) -
            exp(-t_max / td) * (cos(omega * t_max) + static_cast<double>(1.0) / (omega * td) * sin(omega * t_max));

    // Weber number without velocity
    const double We_i_wo_v = rho_a * L / sigma;

    // Equation (8)
    const double y_coeff = (C_F * We_i_wo_v * c_def) / (C_k * C_b);

    const double n_full = 38;
    const double n_full_23 = n_full * 2.0 / 3.0;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)numParticles; i++) {
            const Vector3D &vi = m_model->getVelocity(i);
            Vector3D v_i_rel = va - vi;
            const double vi_rel_square = v_i_rel.lengthSquared();
            const double vi_rel_norm = sqrt(vi_rel_square);

            Vector3D v_i_rel_n = v_i_rel;
            if (vi_rel_norm <= 1.0e-6) continue;
            // Else.
            v_i_rel_n = v_i_rel_n * (1.0 / vi_rel_norm);

            // Equation (8)
            const double y_i_max = std::min(vi_rel_square * y_coeff, static_cast<double>(1.0));

            const double Re_i =
                    static_cast<double>(2.0) * std::max((rho_a * vi_rel_norm * L) / mu_a, static_cast<double>(0.1));

            double C_Di_sphere;
            if (Re_i <= 1000.0)
                C_Di_sphere = static_cast<double>(24.0) / Re_i *
                              (static_cast<double>(1.0) +
                               static_cast<double>(1.0 / 6.0) * pow(Re_i, static_cast<double>(2.0 / 3.0)));
            else
                C_Di_sphere = static_cast<double>(0.424);

            // Equation (9)
            const double C_Di_Liu = C_Di_sphere * (static_cast<double>(1.0) + static_cast<double>(2.632) * y_i_max);

            unsigned int numNeighbors = 0;
            for (unsigned int pid = 0; pid < nFluids; pid++)
                numNeighbors += sim->numberOfNeighbors(fluidModelIndex, pid, i);

            if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
                for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++)
                    numNeighbors += sim->numberOfNeighbors(fluidModelIndex, pid, i);
            }
            // if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
            //{
            //	forall_density_maps(
            //		// ToDo
            //	);
            // }
            else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019) {
                forall_volume_maps(numNeighbors += (unsigned int)(Vj / fluidParticleVolume) + 1;)
            }

            // Equation (10)
            double C_Di;
            const double factor_n = std::min(n_full_23, (double)numNeighbors) / n_full_23;
            if (numNeighbors == 0)
                C_Di = C_Di_Liu;
            else
                C_Di = (static_cast<double>(1.0) - factor_n) * C_Di_Liu + factor_n;

            // Equation (12)
            const double h1 = (L + C_b * L * y_i_max);
            const double A_i_droplet = pi * h1 * h1;

            // Equation (13)
            const double A_i_unoccluded = (static_cast<double>(1.0) - factor_n) * A_i_droplet + factor_n * diam * diam;

            //////////////////////////////////////////////////////////////////////////
            // Fluid
            //////////////////////////////////////////////////////////////////////////
            double max_v_x = 0.0;
            const Vector3D &xi = m_model->getPosition(i);
            forall_fluid_neighbors_in_same_phase(Vector3D xixj = xi - xj; xixj.normalize();
                                                 const double x_v = v_i_rel_n.dot(xixj);
                                                 max_v_x = std::max(max_v_x, x_v);)
                    // Equation (15)
                    const double w_i = std::max(static_cast<double>(0.0),
                                                std::min(static_cast<double>(1.0), static_cast<double>(1.0) - max_v_x));

            // Equation (14)
            const double A_i = w_i * A_i_unoccluded;

            // Drag force. Additionally dividing by mass to get acceleration.
            Vector3D &ai = m_model->getAcceleration(i);
            ai += m_dragCoefficient * static_cast<double>(0.5) / m_model->getMass(i) * rho_a * (v_i_rel * vi_rel_norm) *
                  C_Di * A_i;
        }
    }
}

#endif

void DragForce_Gissler2017::reset() {}

}  // namespace vox::flex