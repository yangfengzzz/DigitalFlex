//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/boundary_model.h"

#include <iostream>

namespace vox::flex {

BoundaryModel::BoundaryModel() : m_forcePerThread(), m_torquePerThread() {}

BoundaryModel::~BoundaryModel() {
    m_forcePerThread.clear();
    m_torquePerThread.clear();

    delete m_rigidBody;
}

void BoundaryModel::reset() {
    for (int j = 0; j < m_forcePerThread.size(); j++) {
        m_forcePerThread[j].setZero();
        m_torquePerThread[j].setZero();
    }
}

void BoundaryModel::getForceAndTorque(Vector3D &force, Vector3D &torque) {
    force.setZero();
    torque.setZero();
    for (int j = 0; j < m_forcePerThread.size(); j++) {
        force += m_forcePerThread[j];
        torque += m_torquePerThread[j];
    }
}

void BoundaryModel::clearForceAndTorque() {
    for (int j = 0; j < m_forcePerThread.size(); j++) {
        m_forcePerThread[j].setZero();
        m_torquePerThread[j].setZero();
    }
}

}  // namespace vox::flex
