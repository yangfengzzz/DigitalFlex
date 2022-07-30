// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.force/zhu_bridson_points_to_implicit2.h"

#include "vox.base/logging.h"
#include "vox.force/fmm_level_set_solver2.h"
#include "vox.geometry/particle_system_data2.h"

using namespace vox;

inline double k(double s) { return std::max(0.0, cubic(1.0 - s * s)); }

ZhuBridsonPointsToImplicit2::ZhuBridsonPointsToImplicit2(double kernelRadius, double cutOffThreshold, bool isOutputSdf)
    : _kernelRadius(kernelRadius), _cutOffThreshold(cutOffThreshold), _isOutputSdf(isOutputSdf) {}

void ZhuBridsonPointsToImplicit2::convert(const ConstArrayAccessor1<Point2D>& points, ScalarGrid2* output) const {
    if (output == nullptr) {
        LOGW("Null scalar grid output pointer provided.")
        return;
    }

    const auto res = output->resolution();
    if (res.x * res.y == 0) {
        LOGW("Empty grid is provided.")
        return;
    }

    const auto bbox = output->boundingBox();
    if (bbox.isEmpty()) {
        LOGW("Empty domain is provided.")
        return;
    }

    ParticleSystemData2 particles;
    particles.addParticles(points);
    particles.buildNeighborSearcher(_kernelRadius);

    const auto neighborSearcher = particles.neighborSearcher();
    const double isoContValue = _cutOffThreshold * _kernelRadius;

    auto temp = output->clone();
    temp->fill([&](const Point2D& x) -> double {
        Point2D xAvg;
        double wSum = 0.0;
        const auto func = [&](size_t, const Point2D& xi) {
            const double wi = k((x - xi).length() / _kernelRadius);

            wSum += wi;
            xAvg += Vector2D(xi.x * wi, xi.y * wi);
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
        FmmLevelSetSolver2 solver;
        solver.reinitialize(*temp, kMaxD, output);
    } else {
        temp->swap(output);
    }
}
