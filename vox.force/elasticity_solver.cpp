//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/elasticity_solver.h"

namespace vox::flex {

int ElasticityBase::YOUNGS_MODULUS = -1;
int ElasticityBase::POISSON_RATIO = -1;
int ElasticityBase::FIXED_BOX_MIN = -1;
int ElasticityBase::FIXED_BOX_MAX = -1;

ElasticityBase::ElasticityBase(FluidModel* model)
    : NonPressureForceBase(model),
      m_youngsModulus(static_cast<double>(100000.0)),
      m_poissonRatio(static_cast<double>(0.3)) {
    m_fixedBoxMin.setZero();
    m_fixedBoxMax.setZero();
}

ElasticityBase::~ElasticityBase() = default;

void ElasticityBase::initParameters() {}

/** Mark all particles in the bounding box as fixed.
 */
void ElasticityBase::determineFixedParticles() {
    const unsigned int numParticles = m_model->numActiveParticles();

    if (!m_fixedBoxMin.isEqual(Vector3D()) || !m_fixedBoxMax.isEqual(Vector3D())) {
        for (int i = 0; i < (int)numParticles; i++) {
            const Vector3D& x = m_model->getPosition(i);
            if ((x.x > m_fixedBoxMin.x) && (x.y > m_fixedBoxMin.y) && (x.z > m_fixedBoxMin.z) &&
                (x.x < m_fixedBoxMax.x) && (x.y < m_fixedBoxMax.y) && (x.z < m_fixedBoxMax.z)) {
                m_model->setParticleState(i, ParticleState::Fixed);
            }
        }
    }
}

}  // namespace vox::flex
