// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "manual_tests.h"
#include "vox.force/particle_system_solver2.h"
#include "vox.force/rigid_body_collider2.h"
#include "vox.geometry/constant_vector_field2.h"
#include "vox.geometry/plane2.h"
#include "vox.geometry/point_particle_emitter2.h"

using namespace vox;

JET_TESTS(ParticleSystemSolver2);

JET_BEGIN_TEST_F(ParticleSystemSolver2, Update) {
    Plane2Ptr plane = std::make_shared<Plane2>(Vector2D(0, 1), Point2D());
    RigidBodyCollider2Ptr collider = std::make_shared<RigidBodyCollider2>(plane);
    ConstantVectorField2Ptr wind = std::make_shared<ConstantVectorField2>(Vector2D(1, 0));

    ParticleSystemSolver2 solver;
    solver.setCollider(collider);
    solver.setWind(wind);

    ParticleSystemData2Ptr particles = solver.particleSystemData();
    PointParticleEmitter2Ptr emitter =
            std::make_shared<PointParticleEmitter2>(Point2D(0, 3), Vector2D(0, 1), 5.0, 45.0);
    emitter->setMaxNumberOfNewParticlesPerSecond(100);
    solver.setEmitter(emitter);

    saveParticleDataXy(particles, 0);

    for (Frame frame(0, 1.0 / 60.0); frame.index < 360; frame.advance()) {
        solver.update(frame);

        saveParticleDataXy(particles, frame.index);
    }
}
JET_END_TEST_F
