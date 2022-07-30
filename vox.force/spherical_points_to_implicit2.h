//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/points_to_implicit2.h"

namespace vox {

//!
//! \brief 2-D points-to-implicit converter based on simple sphere model.
//!
class SphericalPointsToImplicit2 final : public PointsToImplicit2 {
public:
    //! Constructs the converter with given sphere radius.
    explicit SphericalPointsToImplicit2(double radius = 1.0, bool isOutputSdf = true);

    //! Converts the given points to implicit surface scalar field.
    void convert(const ConstArrayAccessor1<Point2D>& points, ScalarGrid2* output) const override;

private:
    double _radius = 1.0;
    bool _isOutputSdf = true;
};

//! Shared pointer type for SphericalPointsToImplicit2.
typedef std::shared_ptr<SphericalPointsToImplicit2> SphericalPointsToImplicit2Ptr;

}  // namespace vox