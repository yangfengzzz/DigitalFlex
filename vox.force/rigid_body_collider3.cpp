// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.force/rigid_body_collider3.h"

using namespace vox;

RigidBodyCollider3::RigidBodyCollider3(const Surface3Ptr& surface) { setSurface(surface); }

RigidBodyCollider3::RigidBodyCollider3(const Surface3Ptr& surface,
                                       const Vector3D& linearVelocity_,
                                       const Vector3D& angularVelocity_)
    : linearVelocity(linearVelocity_), angularVelocity(angularVelocity_) {
    setSurface(surface);
}

Vector3D RigidBodyCollider3::velocityAt(const Point3D& point) const {
    Point3D r = point - surface()->transform.translation();
    return linearVelocity + angularVelocity.cross(Vector3D(r.x, r.y, r.z));
}

RigidBodyCollider3::Builder RigidBodyCollider3::builder() { return {}; }

RigidBodyCollider3::Builder& RigidBodyCollider3::Builder::withSurface(const Surface3Ptr& surface) {
    _surface = surface;
    return *this;
}

RigidBodyCollider3::Builder& RigidBodyCollider3::Builder::withLinearVelocity(const Vector3D& linearVelocity) {
    _linearVelocity = linearVelocity;
    return *this;
}

RigidBodyCollider3::Builder& RigidBodyCollider3::Builder::withAngularVelocity(const Vector3D& angularVelocity) {
    _angularVelocity = angularVelocity;
    return *this;
}

RigidBodyCollider3 RigidBodyCollider3::Builder::build() const {
    return RigidBodyCollider3(_surface, _linearVelocity, _angularVelocity);
}

RigidBodyCollider3Ptr RigidBodyCollider3::Builder::makeShared() const {
    return {new RigidBodyCollider3(_surface, _linearVelocity, _angularVelocity),
            [](RigidBodyCollider3* obj) { delete obj; }};
}
