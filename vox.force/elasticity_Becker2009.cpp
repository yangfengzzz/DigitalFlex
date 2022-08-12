//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/elasticity_Becker2009.h"

#include "vox.force/simulation.h"
#include "vox.geometry/matrix_utils.h"

namespace vox::flex {

int Elasticity_Becker2009::ALPHA = -1;

Elasticity_Becker2009::Elasticity_Becker2009(FluidModel *model) : ElasticityBase(model) {
    const unsigned int numParticles = model->numActiveParticles();
    m_restVolumes.resize(numParticles);
    m_current_to_initial_index.resize(numParticles);
    m_initial_to_current_index.resize(numParticles);
    m_initialNeighbors.resize(numParticles);
    m_rotations.resize(numParticles, Matrix3x3D::makeIdentity());
    m_stress.resize(numParticles);
    m_F.resize(numParticles);
    m_alpha = 0.0;

    initValues();

    model->addField({"rest volume", FieldType::Scalar,
                     [&](const unsigned int i) -> Real * { return &m_restVolumes[i]; }, true});
    model->addField(
            {"rotation", FieldType::Matrix3, [&](const unsigned int i) -> Real * { return &m_rotations[i](0, 0); }});
    model->addField({"stress", FieldType::Vector6, [&](const unsigned int i) -> Real * { return &m_stress[i][0]; }});
    model->addField({"deformation gradient", FieldType::Matrix3,
                     [&](const unsigned int i) -> Real * { return &m_F[i](0, 0); }});
}

Elasticity_Becker2009::~Elasticity_Becker2009() {
    m_model->removeFieldByName("rest volume");
    m_model->removeFieldByName("rotation");
    m_model->removeFieldByName("stress");
    m_model->removeFieldByName("deformation gradient");
}

void Elasticity_Becker2009::initParameters() {
    ElasticityBase::initParameters();

    ALPHA = createNumericParameter("alpha", "Zero-energy modes suppression", &m_alpha);
    setGroup(ALPHA, "Elasticity");
    setDescription(ALPHA, "Coefficient for zero-energy modes suppression method");
    auto *rparam = static_cast<utility::DoubleParameter *>(getParameter(ALPHA));
    rparam->setMinValue(0.0);
}

void Elasticity_Becker2009::initValues() {
    Simulation *sim = Simulation::getCurrent();
    sim->getNeighborhoodSearch()->find_neighbors();

    FluidModel *model = m_model;
    const unsigned int numParticles = model->numActiveParticles();
    const unsigned int fluidModelIndex = model->getPointSetIndex();

// Store the neighbors in the reference configurations and
// compute the volume of each particle in rest state
#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)numParticles; i++) {
            m_current_to_initial_index[i] = i;
            m_initial_to_current_index[i] = i;

            // only neighbors in same phase will influence elasticity
            const unsigned int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);
            m_initialNeighbors[i].resize(numNeighbors);
            for (unsigned int j = 0; j < numNeighbors; j++)
                m_initialNeighbors[i][j] = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);

            // compute volume
            Real density = model->getMass(i) * sim->W_zero();
            const Vector3D &xi = model->getPosition(i);
            forall_fluid_neighbors_in_same_phase(density += model->getMass(neighborIndex) * sim->W(xi - xj);)
                    m_restVolumes[i] = model->getMass(i) / density;
        }
    }

    // mark all particles in the bounding box as fixed
    determineFixedParticles();
}

void Elasticity_Becker2009::step() {
    computeRotations();
    computeStress();
    computeForces();
}

void Elasticity_Becker2009::reset() { initValues(); }

void Elasticity_Becker2009::performNeighborhoodSearchSort() {
    const unsigned int numPart = m_model->numActiveParticles();
    if (numPart == 0) return;

    Simulation *sim = Simulation::getCurrent();
    auto const &d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
    d.sort_field(&m_restVolumes[0]);
    d.sort_field(&m_current_to_initial_index[0]);

    for (unsigned int i = 0; i < numPart; i++) m_initial_to_current_index[m_current_to_initial_index[i]] = i;
}

void Elasticity_Becker2009::computeRotations() {
    Simulation *sim = Simulation::getCurrent();
    const unsigned int numParticles = m_model->numActiveParticles();
    FluidModel *model = m_model;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)numParticles; i++) {
            const unsigned int i0 = m_current_to_initial_index[i];
            const Vector3D &xi = m_model->getPosition(i);
            const Vector3D &xi0 = m_model->getPosition0(i0);
            Matrix3x3D Apq;

            const size_t numNeighbors = m_initialNeighbors[i0].size();

            //////////////////////////////////////////////////////////////////////////
            // Fluid
            //////////////////////////////////////////////////////////////////////////
            for (unsigned int j = 0; j < numNeighbors; j++) {
                const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
                // get initial neighbor index considering the current particle order
                const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

                const Vector3D &xj = model->getPosition(neighborIndex);
                const Vector3D &xj0 = m_model->getPosition0(neighborIndex0);
                const Vector3D xj_xi = xj - xi;
                const Vector3D xj_xi_0 = xj0 - xi0;
                Apq += m_model->getMass(neighborIndex) * sim->W(xj_xi_0) * Matrix3x3D::makeTensorMatrix(xj_xi, xj_xi_0);
            }

            // Vector3D sigma;
            // Matrix3x3D U, VT;
            // svdWithInversionHandling(Apq, sigma, U, VT);
            // m_rotations[i] = U * VT;
            QuaternionD q(m_rotations[i]);
            extractRotation(Apq, q, 10);
            m_rotations[i] = q.matrix3();
        }
    }
}

void Elasticity_Becker2009::computeStress() {
    Simulation *sim = Simulation::getCurrent();
    const unsigned int numParticles = m_model->numActiveParticles();
    FluidModel *model = m_model;

    // Elasticity tensor
    Matrix<double, 6, 6> C;
    const Real factor = m_youngsModulus / ((static_cast<Real>(1.0) + m_poissonRatio) *
                                           (static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio));
    C(0, 0) = C(1, 1) = C(2, 2) = factor * (static_cast<Real>(1.0) - m_poissonRatio);
    C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = factor * (m_poissonRatio);
    C(3, 3) = C(4, 4) = C(5, 5) =
            factor * static_cast<Real>(0.5) * (static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio);

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)numParticles; i++) {
            if (model->getParticleState(i) == ParticleState::Active) {
                const unsigned int i0 = m_current_to_initial_index[i];
                const Vector3D &xi = m_model->getPosition(i);
                const Vector3D &xi0 = m_model->getPosition0(i0);

                Matrix3x3D nablaU;
                const size_t numNeighbors = m_initialNeighbors[i0].size();

                //////////////////////////////////////////////////////////////////////////
                // Fluid
                //////////////////////////////////////////////////////////////////////////
                for (unsigned int j = 0; j < numNeighbors; j++) {
                    const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
                    // get initial neighbor index considering the current particle order
                    const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

                    const Vector3D &xj = model->getPosition(neighborIndex);
                    const Vector3D &xj0 = m_model->getPosition0(neighborIndex0);

                    const Vector3D xj_xi = xj - xi;
                    const Vector3D xj_xi_0 = xj0 - xi0;

                    const Vector3D uji = m_rotations[i].transposed() * xj_xi - xj_xi_0;
                    // subtract because kernel gradient is taken in direction of xji0 instead of xij0
                    nablaU -= Matrix3x3D::makeTensorMatrix(m_restVolumes[neighborIndex] * uji, sim->gradW(xj_xi_0));
                }
                m_F[i] = nablaU + Matrix3x3D::makeIdentity();

                // compute Cauchy strain: epsilon = 0.5 (nabla u + nabla u^T)
                Vector<double, 6> strain;
                strain[0] = nablaU(0, 0);                                            // \epsilon_{00}
                strain[1] = nablaU(1, 1);                                            // \epsilon_{11}
                strain[2] = nablaU(2, 2);                                            // \epsilon_{22}
                strain[3] = static_cast<Real>(0.5) * (nablaU(0, 1) + nablaU(1, 0));  // \epsilon_{01}
                strain[4] = static_cast<Real>(0.5) * (nablaU(0, 2) + nablaU(2, 0));  // \epsilon_{02}
                strain[5] = static_cast<Real>(0.5) * (nablaU(1, 2) + nablaU(2, 1));  // \epsilon_{12}

                // stress = C * epsilon
                m_stress[i] = C * strain;
            } else
                m_stress[i].setZero();
        }
    }
}

void Elasticity_Becker2009::computeForces() {
    Simulation *sim = Simulation::getCurrent();
    const unsigned int numParticles = m_model->numActiveParticles();
    FluidModel *model = m_model;

#pragma omp parallel default(shared)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < (int)numParticles; i++) {
            if (model->getParticleState(i) == ParticleState::Active) {
                const unsigned int i0 = m_current_to_initial_index[i];
                const Vector3D &xi0 = m_model->getPosition0(i0);

                const size_t numNeighbors = m_initialNeighbors[i0].size();
                Vector3D fi;
                fi.setZero();

                //////////////////////////////////////////////////////////////////////////
                // Fluid
                //////////////////////////////////////////////////////////////////////////
                for (unsigned int j = 0; j < numNeighbors; j++) {
                    const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
                    // get initial neighbor index considering the current particle order
                    const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

                    const Vector3D &xj0 = m_model->getPosition0(neighborIndex0);

                    const Vector3D xj_xi_0 = xj0 - xi0;
                    const Vector3D gradW0 = sim->gradW(xj_xi_0);

                    const Vector3D dji = m_restVolumes[i] * gradW0;
                    const Vector3D dij = -m_restVolumes[neighborIndex] * gradW0;

                    Vector3D sdji, sdij;
                    symMatTimesVec(m_stress[neighborIndex], dji, sdji);
                    symMatTimesVec(m_stress[i], dij, sdij);

                    const Vector3D fij = -m_restVolumes[neighborIndex] * sdji;
                    const Vector3D fji = -m_restVolumes[i] * sdij;

                    fi += m_rotations[neighborIndex] * fij - m_rotations[i] * fji;
                }
                fi = 0.5 * fi;

                if (m_alpha != 0.0) {
                    //////////////////////////////////////////////////////////////////////////
                    // Ganzenmüller, G.C. 2015. An hourglass control algorithm for Lagrangian
                    // Smooth Particle Hydrodynamics. Computer Methods in Applied Mechanics and
                    // Engineering 286, 87.106.
                    //////////////////////////////////////////////////////////////////////////
                    Vector3D fi_hg;
                    fi_hg.setZero();
                    const Vector3D &xi = m_model->getPosition(i);
                    for (unsigned int j = 0; j < numNeighbors; j++) {
                        const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
                        // get initial neighbor index considering the current particle order
                        const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

                        const Vector3D &xj = model->getPosition(neighborIndex);
                        const Vector3D &xj0 = m_model->getPosition0(neighborIndex0);

                        // Note: Ganzenm�ller defines xij = xj-xi
                        const Vector3D xi_xj = -(xi - xj);
                        const Real xixj_l = xi_xj.length();
                        if (xixj_l > 1.0e-6) {
                            // Note: Ganzenm�ller defines xij = xj-xi
                            const Vector3D xi_xj_0 = -(xi0 - xj0);
                            const Real xixj0_l2 = xi_xj_0.lengthSquared();
                            const Real W0 = sim->W(xi_xj_0);

                            const Vector3D xij_i = m_F[i] * m_rotations[i] * xi_xj_0;
                            const Vector3D xji_j = -m_F[neighborIndex] * m_rotations[neighborIndex] * xi_xj_0;
                            const Vector3D epsilon_ij_i = xij_i - xi_xj;
                            const Vector3D epsilon_ji_j = xji_j + xi_xj;

                            const Real delta_ij_i = epsilon_ij_i.dot(xi_xj) / xixj_l;
                            const Real delta_ji_j = -epsilon_ji_j.dot(xi_xj) / xixj_l;

                            fi_hg -= m_restVolumes[neighborIndex] * W0 / xixj0_l2 * (delta_ij_i + delta_ji_j) * xi_xj /
                                     xixj_l;
                        }
                    }
                    fi_hg *= m_alpha * m_youngsModulus * m_restVolumes[i];
                    model->getAcceleration(i) += fi_hg / model->getMass(i);
                }

                // elastic acceleration
                Vector3D &ai = model->getAcceleration(i);
                ai += fi / model->getMass(i);
            }
        }
    }
}

void Elasticity_Becker2009::saveState(BinaryFileWriter &binWriter) {
    binWriter.writeBuffer((char *)m_current_to_initial_index.data(),
                          m_current_to_initial_index.size() * sizeof(unsigned int));
    binWriter.writeBuffer((char *)m_initial_to_current_index.data(),
                          m_initial_to_current_index.size() * sizeof(unsigned int));
}

void Elasticity_Becker2009::loadState(BinaryFileReader &binReader) {
    binReader.readBuffer((char *)m_current_to_initial_index.data(),
                         m_current_to_initial_index.size() * sizeof(unsigned int));
    binReader.readBuffer((char *)m_initial_to_current_index.data(),
                         m_initial_to_current_index.size() * sizeof(unsigned int));
}

}  // namespace vox::flex