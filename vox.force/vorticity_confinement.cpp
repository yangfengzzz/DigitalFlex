//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/vorticity_confinement.h"

#include "vox.force/simulation.h"

namespace vox::flex {

VorticityConfinement::VorticityConfinement(FluidModel *model) : VorticityBase(model) {
    m_omega.resize(model->numParticles(), Vector3D());
    m_normOmega.resize(model->numParticles(), 0.0);

    model->addField(
            {"angular velocity", FieldType::Vector3, [&](const unsigned int i) -> Real * { return &m_omega[i][0]; }});
}

VorticityConfinement::~VorticityConfinement() {
    m_model->removeFieldByName("angular velocity");

    m_omega.clear();
    m_normOmega.clear();
}

void VorticityConfinement::step() {
    Simulation *sim = Simulation::getCurrent();
    const unsigned int numParticles = m_model->numActiveParticles();
    const unsigned int fluidModelIndex = m_model->getPointSetIndex();
    FluidModel *model = m_model;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)numParticles; i++) {
            const Vector3D &xi = m_model->getPosition(i);
            const Vector3D &vi = m_model->getVelocity(i);
            Vector3D &omegai = m_omega[i];
            omegai.setZero();
            const Real density_i = m_model->getDensity(i);

            //////////////////////////////////////////////////////////////////////////
            // Fluid
            //////////////////////////////////////////////////////////////////////////
            forall_fluid_neighbors_in_same_phase(const Vector3D &vj = m_model->getVelocity(neighborIndex);
                                                 const Real density_j = m_model->getDensity(neighborIndex);
                                                 const Real density_j2 = density_j * density_j;
                                                 const Vector3D gradW = sim->gradW(xi - xj);

                                                 omegai -= m_model->getMass(neighborIndex) / density_i *
                                                           (vi - vj).cross(gradW);) Real &normOmegai = m_normOmega[i];
            normOmegai = omegai.length();
        }
    }

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)numParticles; i++) {
            const Vector3D &xi = m_model->getPosition(i);
            Vector3D &ai = m_model->getAcceleration(i);
            const Real density_i = m_model->getDensity(i);
            const Vector3D &omegai = m_omega[i];

            Vector3D etai;
            etai.setZero();
            //////////////////////////////////////////////////////////////////////////
            // Fluid
            //////////////////////////////////////////////////////////////////////////
            forall_fluid_neighbors_in_same_phase(
                    const Real density_j = m_model->getDensity(neighborIndex);
                    const Vector3D gradW = sim->gradW(xi - xj); Real &normOmegaj = m_normOmega[neighborIndex];
                    etai += m_model->getMass(neighborIndex) / density_i * normOmegaj * gradW;)

                    etai.normalize();
            ai += m_vorticityCoeff * etai.cross(omegai);
        }
    }
}

void VorticityConfinement::reset() {}

void VorticityConfinement::performNeighborhoodSearchSort() {}

}  // namespace vox::flex