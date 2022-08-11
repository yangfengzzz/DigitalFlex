//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/boundary_simulator.h"

#include "vox.force/boundary_model.h"
#include "vox.force/simulation.h"

namespace vox::flex {

void BoundarySimulator::updateBoundaryForces() {
    Simulation *sim = Simulation::getCurrent();
    const unsigned int nObjects = sim->numberOfBoundaryModels();
    for (unsigned int i = 0; i < nObjects; i++) {
        BoundaryModel *bm = sim->getBoundaryModel(i);
        RigidBodyObject *rbo = bm->getRigidBodyObject();
        if (rbo->isDynamic()) {
            Vector3D force, torque;
            bm->getForceAndTorque(force, torque);
            rbo->addForce(force);
            rbo->addTorque(torque);
            bm->clearForceAndTorque();
        }
    }
}
}  // namespace vox::flex