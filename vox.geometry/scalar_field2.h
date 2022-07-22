// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_JET_SCALAR_FIELD2_H_
#define INCLUDE_JET_SCALAR_FIELD2_H_

#include <functional>
#include <memory>

#include "vox.geometry/field2.h"
#include "vox.math/point2.h"

namespace vox {

//! Abstract base class for 2-D scalar field.
class ScalarField2 : public Field2 {
public:
    //! Default constructor.
    ScalarField2();

    //! Default destructor.
    virtual ~ScalarField2();

    //! Returns sampled value at given position \p x.
    virtual double sample(const Point2D &x) const = 0;

    //! Returns gradient vector at given position \p x.
    virtual Vector2D gradient(const Point2D &x) const;

    //! Returns Laplacian at given position \p x.
    virtual double laplacian(const Point2D &x) const;

    //! Returns sampler function object.
    virtual std::function<double(const Point2D &)> sampler() const;
};

//! Shared pointer for the ScalarField2 type.
typedef std::shared_ptr<ScalarField2> ScalarField2Ptr;

}  // namespace vox

#endif  // INCLUDE_JET_SCALAR_FIELD2_H_
