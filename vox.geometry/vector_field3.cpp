// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.geometry/vector_field3.h"

using namespace vox;

VectorField3::VectorField3() {}

VectorField3::~VectorField3() {}

double VectorField3::divergence(const Point3D &) const { return 0.0; }

Vector3D VectorField3::curl(const Point3D &) const { return Vector3D(); }

std::function<Vector3D(const Point3D &)> VectorField3::sampler() const {
    const VectorField3 *self = this;
    return [self](const Point3D &x) -> Vector3D { return self->sample(x); };
}
