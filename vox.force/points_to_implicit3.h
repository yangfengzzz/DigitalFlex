//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <memory>

#include "vox.geometry/array_accessor1.h"
#include "vox.geometry/scalar_grid3.h"
#include "vox.math/point3.h"

namespace vox {

//! Abstract base class for 3-D points-to-implicit converters.
class PointsToImplicit3 {
public:
    //! Default constructor.
    PointsToImplicit3();

    //! Default destructor.
    virtual ~PointsToImplicit3();

    //! Converts the given points to implicit surface scalar field.
    virtual void convert(const ConstArrayAccessor1<Point3D>& points, ScalarGrid3* output) const = 0;
};

//! Shared pointer for the PointsToImplicit3 type.
typedef std::shared_ptr<PointsToImplicit3> PointsToImplicit3Ptr;

}  // namespace vox