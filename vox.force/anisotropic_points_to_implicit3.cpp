// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.force/anisotropic_points_to_implicit3.h"

#include "vox.base/logging.h"
#include "vox.force/fmm_level_set_solver3.h"
#include "vox.force/sph_kernels3.h"
#include "vox.force/sph_system_data3.h"
#include "vox.geometry/point_kdtree_searcher3.h"
#include "vox.geometry/svd.h"

using namespace vox;

inline double p(double distance) {
    const double distanceSquared = distance * distance;

    if (distanceSquared >= 1.0) {
        return 0.0;
    } else {
        const double x = 1.0 - distanceSquared;
        return x * x * x;
    }
}

inline double wij(double distance, double r) {
    if (distance < r) {
        return 1.0 - cubic(distance / r);
    } else {
        return 0.0;
    }
}

inline Matrix3x3D vvt(const Vector3D& v) {
    return {v.x * v.x, v.x * v.y, v.x * v.z, v.y * v.x, v.y * v.y, v.y * v.z, v.z * v.x, v.z * v.y, v.z * v.z};
}

inline double w(const Vector3D& r, const Matrix3x3D& g, double gDet) {
    static const double sigma = 315.0 / (64 * kPiD);
    return sigma * gDet * p((g * r).length());
}

//

AnisotropicPointsToImplicit3::AnisotropicPointsToImplicit3(double kernelRadius,
                                                           double cutOffDensity,
                                                           double positionSmoothingFactor,
                                                           size_t minNumNeighbors,
                                                           bool isOutputSdf)
    : _kernelRadius(kernelRadius),
      _cutOffDensity(cutOffDensity),
      _positionSmoothingFactor(positionSmoothingFactor),
      _minNumNeighbors(minNumNeighbors),
      _isOutputSdf(isOutputSdf) {}

void AnisotropicPointsToImplicit3::convert(const ConstArrayAccessor1<Point3D>& points, ScalarGrid3* output) const {
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

    LOGI("Start converting points to implicit surface.")

    const double h = _kernelRadius;
    const double invH = 1 / h;
    const double r = 2.0 * h;

    // Mean estimator for cov. mat.
    const auto meanNeighborSearcher = PointKdTreeSearcher3::builder().makeShared();
    meanNeighborSearcher->build(points);

    LOGI("Built neighbor searcher.")

    SphSystemData3 meanParticles;
    meanParticles.addParticles(points);
    meanParticles.setNeighborSearcher(meanNeighborSearcher);
    meanParticles.setKernelRadius(r);

    // Compute G and xMean
    std::vector<Matrix3x3D> gs(points.size());
    Array1<Point3D> xMeans(points.size());

    parallelFor(kZeroSize, points.size(), [&](size_t i) {
        const auto& x = points[i];

        // Compute xMean
        Point3D xMean;
        double wSum = 0.0;
        size_t numNeighbors = 0;
        const auto getXMean = [&](size_t, const Point3D& xj) {
            const double wj = wij((x - xj).length(), r);
            wSum += wj;
            xMean += Vector3D(wj * xj.x, wj * xj.y, wj * xj.z);
            ++numNeighbors;
        };
        meanNeighborSearcher->forEachNearbyPoint(x, r, getXMean);

        VOX_ASSERT(wSum > 0.0);
        xMean /= wSum;

        xMeans[i] = Point3D(lerp(x.x, xMean.x, _positionSmoothingFactor), lerp(x.y, xMean.y, _positionSmoothingFactor),
                            lerp(x.z, xMean.z, _positionSmoothingFactor));

        if (numNeighbors < _minNumNeighbors) {
            const auto g = Matrix3x3D::makeScaleMatrix(invH, invH, invH);
            gs[i] = g;
        } else {
            // Compute covariance matrix
            // We start with small scale matrix (h*h) in order to
            // prevent zero covariance matrix when points are all
            // perfectly lined up.
            auto cov = Matrix3x3D::makeScaleMatrix(h * h, h * h, h * h);
            wSum = 0.0;
            const auto getCov = [&](size_t, const Point3D& xj) {
                const double wj = wij((xMean - xj).length(), r);
                wSum += wj;
                cov += wj * vvt(xj - xMean);
            };
            meanNeighborSearcher->forEachNearbyPoint(x, r, getCov);

            cov /= wSum;

            // SVD
            Matrix3x3D u;
            Vector3D v;
            Matrix3x3D w;
            svd(cov, u, v, w);

            // Take off the sign
            v.x = std::fabs(v.x);
            v.y = std::fabs(v.y);
            v.z = std::fabs(v.z);

            // Constrain Sigma
            const double maxSingularVal = v.max();
            const double kr = 4.0;
            v.x = std::max(v.x, maxSingularVal / kr);
            v.y = std::max(v.y, maxSingularVal / kr);
            v.z = std::max(v.z, maxSingularVal / kr);

            const auto invSigma = Matrix3x3D::makeScaleMatrix(1.0 / v);

            // Compute G
            const double scale = std::pow(v.x * v.y * v.z, 1.0 / 3.0);  // volume preservation
            const Matrix3x3D g = invH * scale * (w * invSigma * u.transposed());
            gs[i] = g;
        }
    });

    LOGI("Computed G and means.")

    // SPH estimator
    meanParticles.setKernelRadius(h);
    meanParticles.updateDensities();
    const auto d = meanParticles.densities();
    const double m = meanParticles.mass();

    PointKdTreeSearcher3 meanNeighborSearcher2;
    meanNeighborSearcher2.build(xMeans.constAccessor());

    // Compute SDF
    auto temp = output->clone();
    temp->fill([&](const Point3D& x) {
        double sum = 0.0;
        meanNeighborSearcher2.forEachNearbyPoint(x, r, [&](size_t i, const Point3D& neighborPosition) {
            sum += m / d[i] * w(neighborPosition - x, gs[i], gs[i].determinant());
        });

        return _cutOffDensity - sum;
    });

    LOGI("Computed SDF.")

    if (_isOutputSdf) {
        FmmLevelSetSolver3 solver;
        solver.reinitialize(*temp, kMaxD, output);

        LOGI("Completed initialization.")
    } else {
        temp->swap(output);
    }

    LOGI("Done converting points to implicit surface.")
}
