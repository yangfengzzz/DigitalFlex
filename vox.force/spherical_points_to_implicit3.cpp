// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.force/spherical_points_to_implicit3.h"

#include "vox.base/logging.h"
#include "vox.force/fmm_level_set_solver3.h"
#include "vox.geometry/particle_system_data3.h"

using namespace vox;

SphericalPointsToImplicit3::SphericalPointsToImplicit3(double radius, bool isOutputSdf)
    : _radius(radius), _isOutputSdf(isOutputSdf) {}

void SphericalPointsToImplicit3::convert(const ConstArrayAccessor1<Point3D>& points, ScalarGrid3* output) const {
    if (output == nullptr) {
        LOGW("Null scalar grid output pointer provided.")
        return;
    }

    const auto res = output->resolution();
    if (res.x * res.y * res.z == 0) {
        LOGW("Empty grid is provided.")
        return;
    }

    const auto bbox = output->boundingBox();
    if (bbox.isEmpty()) {
        LOGW("Empty domain is provided.")
        return;
    }

    ParticleSystemData3 particles;
    particles.addParticles(points);
    particles.buildNeighborSearcher(2.0 * _radius);

    const auto neighborSearcher = particles.neighborSearcher();

    auto temp = output->clone();
    temp->fill([&](const Point3D& x) {
        double minDist = 2.0 * _radius;
        neighborSearcher->forEachNearbyPoint(
                x, 2.0 * _radius, [&](size_t, const Point3D& xj) { minDist = std::min(minDist, (x - xj).length()); });

        return minDist - _radius;
    });

    if (_isOutputSdf) {
        FmmLevelSetSolver3 solver;
        solver.reinitialize(*temp, kMaxD, output);
    } else {
        temp->swap(output);
    }
}
