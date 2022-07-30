//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <memory>

#include "vox.geometry/array_accessor1.h"
#include "vox.geometry/scalar_grid2.h"
#include "vox.math/point2.h"

namespace vox {

//! Abstract base class for 2-D points-to-implicit converters.
class PointsToImplicit2 {
public:
    //! Default constructor.
    PointsToImplicit2();

    //! Default destructor.
    virtual ~PointsToImplicit2();

    //! Converts the given points to implicit surface scalar field.
    virtual void convert(const ConstArrayAccessor1<Point2D>& points, ScalarGrid2* output) const = 0;
};

//! Shared pointer for the PointsToImplicit2 type.
typedef std::shared_ptr<PointsToImplicit2> PointsToImplicit2Ptr;

}  // namespace vox