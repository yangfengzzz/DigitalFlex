// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_JET_VECTOR_FIELD3_H_
#define INCLUDE_JET_VECTOR_FIELD3_H_

#include <functional>
#include <memory>

#include "vox.geometry/field3.h"
#include "vox.math/point3.h"

namespace vox {

//! Abstract base class for 3-D vector field.
class VectorField3 : public Field3 {
public:
    //! Default constructor.
    VectorField3();

    //! Default destructor.
    virtual ~VectorField3();

    //! Returns sampled value at given position \p x.
    virtual Vector3D sample(const Point3D &x) const = 0;

    //! Returns divergence at given position \p x.
    virtual double divergence(const Point3D &x) const;

    //! Returns curl at given position \p x.
    virtual Vector3D curl(const Point3D &x) const;

    //! Returns sampler function object.
    virtual std::function<Vector3D(const Point3D &)> sampler() const;
};

//! Shared pointer for the VectorField3 type.
typedef std::shared_ptr<VectorField3> VectorField3Ptr;

}  // namespace vox

#endif  // INCLUDE_JET_VECTOR_FIELD3_H_
