// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.force/zhu_bridson_points_to_implicit3.h"

#include "vox.base/logging.h"
#include "vox.force/fmm_level_set_solver3.h"
#include "vox.geometry/particle_system_data3.h"

using namespace vox;

inline double k(double s) { return std::max(0.0, cubic(1.0 - s * s)); }

ZhuBridsonPointsToImplicit3::ZhuBridsonPointsToImplicit3(double kernelRadius, double cutOffThreshold, bool isOutputSdf)
    : _kernelRadius(kernelRadius), _cutOffThreshold(cutOffThreshold), _isOutputSdf(isOutputSdf) {}

void ZhuBridsonPointsToImplicit3::convert(const ConstArrayAccessor1<Point3D>& points, ScalarGrid3* output) const {
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
    particles.buildNeighborSearcher(_kernelRadius);

    const auto neighborSearcher = particles.neighborSearcher();
    const double isoContValue = _cutOffThreshold * _kernelRadius;

    auto temp = output->clone();
    temp->fill([&](const Point3D& x) -> double {
        Point3D xAvg;
        double wSum = 0.0;
        const auto func = [&](size_t, const Point3D& xi) {
            const double wi = k((x - xi).length() / _kernelRadius);
            wSum += wi;
            xAvg += Vector3D(wi * xi.x, wi * xi.y, wi * xi.z);
        };
        neighborSearcher->forEachNearbyPoint(x, _kernelRadius, func);

        if (wSum > 0.0) {
            xAvg /= wSum;
            return (x - xAvg).length() - isoContValue;
        } else {
            return output->boundingBox().diagonalLength();
        }
    });

    if (_isOutputSdf) {
        FmmLevelSetSolver3 solver;
        solver.reinitialize(*temp, kMaxD, output);
    } else {
        temp->swap(output);
    }
}
