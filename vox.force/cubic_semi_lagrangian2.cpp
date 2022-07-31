// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.force/cubic_semi_lagrangian2.h"

#include "vox.geometry/array_samplers2.h"

using namespace vox;

CubicSemiLagrangian2::CubicSemiLagrangian2() = default;

std::function<double(const Point2D&)> CubicSemiLagrangian2::getScalarSamplerFunc(const ScalarGrid2& source) const {
    auto sourceSampler =
            CubicArraySampler2<double, double>(source.constDataAccessor(), source.gridSpacing(), source.dataOrigin());
    return sourceSampler.functor();
}

std::function<Vector2D(const Point2D&)> CubicSemiLagrangian2::getVectorSamplerFunc(
        const CollocatedVectorGrid2& source) const {
    auto sourceSampler =
            CubicArraySampler2<Vector2D, double>(source.constDataAccessor(), source.gridSpacing(), source.dataOrigin());
    return sourceSampler.functor();
}

std::function<Vector2D(const Point2D&)> CubicSemiLagrangian2::getVectorSamplerFunc(
        const FaceCenteredGrid2& source) const {
    auto uSourceSampler =
            CubicArraySampler2<double, double>(source.uConstAccessor(), source.gridSpacing(), source.uOrigin());
    auto vSourceSampler =
            CubicArraySampler2<double, double>(source.vConstAccessor(), source.gridSpacing(), source.vOrigin());
    return [uSourceSampler, vSourceSampler](const Point2D& x) {
        return Vector2D(uSourceSampler(x), vSourceSampler(x));
    };
}
