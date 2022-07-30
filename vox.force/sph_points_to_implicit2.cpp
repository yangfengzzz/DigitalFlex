// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.force/sph_points_to_implicit2.h"

#include "vox.base/logging.h"
#include "vox.force/fmm_level_set_solver2.h"
#include "vox.force/sph_system_data2.h"

using namespace vox;

SphPointsToImplicit2::SphPointsToImplicit2(double kernelRadius, double cutOffDensity, bool isOutputSdf)
    : _kernelRadius(kernelRadius), _cutOffDensity(cutOffDensity), _isOutputSdf(isOutputSdf) {}

void SphPointsToImplicit2::convert(const ConstArrayAccessor1<Point2D>& points, ScalarGrid2* output) const {
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

    SphSystemData2 sphParticles;
    sphParticles.addParticles(points);
    sphParticles.setKernelRadius(_kernelRadius);
    sphParticles.buildNeighborSearcher();
    sphParticles.updateDensities();

    Array1<double> constData(sphParticles.numberOfParticles(), 1.0);
    auto temp = output->clone();
    temp->fill([&](const Point2D& x) {
        double d = sphParticles.interpolate(x, constData.constAccessor());
        return _cutOffDensity - d;
    });

    if (_isOutputSdf) {
        FmmLevelSetSolver2 solver;
        solver.reinitialize(*temp, kMaxD, output);
    } else {
        temp->swap(output);
    }
}
