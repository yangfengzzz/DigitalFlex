// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_JET_IMPLICIT_SURFACE2_H_
#define INCLUDE_JET_IMPLICIT_SURFACE2_H_

#include "vox.geometry/surface2.h"

namespace vox {

//! Abstract base class for 2-D implicit surface.
class ImplicitSurface2 : public Surface2 {
public:
    //! Default constructor.
    ImplicitSurface2(const Transform2D &transform = Transform2D(), bool isNormalFlipped = false);

    //! Copy constructor.
    ImplicitSurface2(const ImplicitSurface2 &other);

    //! Default destructor.
    virtual ~ImplicitSurface2();

    //! Returns signed distance from the given point \p otherPoint.
    double signedDistance(const Point2D &otherPoint) const;

protected:
    //! Returns signed distance from the given point \p otherPoint in local
    //! space.
    virtual double signedDistanceLocal(const Point2D &otherPoint) const = 0;

private:
    double closestDistanceLocal(const Point2D &otherPoint) const override;

    bool isInsideLocal(const Point2D &otherPoint) const override;
};

//! Shared pointer type for the ImplicitSurface2.
typedef std::shared_ptr<ImplicitSurface2> ImplicitSurface2Ptr;

}  // namespace vox

#endif  // INCLUDE_JET_IMPLICIT_SURFACE2_H_
