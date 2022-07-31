// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.force/cubic_semi_lagrangian3.h"

#include "vox.geometry/array_samplers3.h"

using namespace vox;

CubicSemiLagrangian3::CubicSemiLagrangian3() = default;

std::function<double(const Point3D&)> CubicSemiLagrangian3::getScalarSamplerFunc(const ScalarGrid3& source) const {
    auto sourceSampler =
            CubicArraySampler3<double, double>(source.constDataAccessor(), source.gridSpacing(), source.dataOrigin());
    return sourceSampler.functor();
}

std::function<Vector3D(const Point3D&)> CubicSemiLagrangian3::getVectorSamplerFunc(
        const CollocatedVectorGrid3& source) const {
    auto sourceSampler =
            CubicArraySampler3<Vector3D, double>(source.constDataAccessor(), source.gridSpacing(), source.dataOrigin());
    return sourceSampler.functor();
}

std::function<Vector3D(const Point3D&)> CubicSemiLagrangian3::getVectorSamplerFunc(
        const FaceCenteredGrid3& source) const {
    auto uSourceSampler =
            CubicArraySampler3<double, double>(source.uConstAccessor(), source.gridSpacing(), source.uOrigin());
    auto vSourceSampler =
            CubicArraySampler3<double, double>(source.vConstAccessor(), source.gridSpacing(), source.vOrigin());
    auto wSourceSampler =
            CubicArraySampler3<double, double>(source.wConstAccessor(), source.gridSpacing(), source.wOrigin());
    return [uSourceSampler, vSourceSampler, wSourceSampler](const Point3D& x) {
        return Vector3D(uSourceSampler(x), vSourceSampler(x), wSourceSampler(x));
    };
}
