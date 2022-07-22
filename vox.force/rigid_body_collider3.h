//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/collider3.h"
#include "vox.math/quaternion.h"

namespace vox {

//!
//! \brief 3-D rigid body collider class.
//!
//! This class implements 3-D rigid body collider. The collider can only take
//! rigid body motion with linear and rotational velocities.
//!
class RigidBodyCollider3 final : public Collider3 {
public:
    class Builder;

    //! Linear velocity of the rigid body.
    Vector3D linearVelocity;

    //! Angular velocity of the rigid body.
    Vector3D angularVelocity;

    //! Constructs a collider with a surface.
    explicit RigidBodyCollider3(const Surface3Ptr& surface);

    //! Constructs a collider with a surface and other parameters.
    RigidBodyCollider3(const Surface3Ptr& surface, const Vector3D& linearVelocity, const Vector3D& angularVelocity);

    //! Returns the velocity of the collider at given \p point.
    [[nodiscard]] Vector3D velocityAt(const Point3D& point) const override;

    //! Returns builder fox RigidBodyCollider3.
    static Builder builder();
};

//! Shared pointer for the RigidBodyCollider3 type.
typedef std::shared_ptr<RigidBodyCollider3> RigidBodyCollider3Ptr;

//!
//! \brief Front-end to create RigidBodyCollider3 objects step by step.
//!
class RigidBodyCollider3::Builder final {
public:
    //! Returns builder with surface.
    Builder& withSurface(const Surface3Ptr& surface);

    //! Returns builder with linear velocity.
    Builder& withLinearVelocity(const Vector3D& linearVelocity);

    //! Returns builder with angular velocity.
    Builder& withAngularVelocity(const Vector3D& angularVelocity);

    //! Builds RigidBodyCollider3.
    [[nodiscard]] RigidBodyCollider3 build() const;

    //! Builds shared pointer of RigidBodyCollider3 instance.
    [[nodiscard]] RigidBodyCollider3Ptr makeShared() const;

private:
    Surface3Ptr _surface;
    Vector3D _linearVelocity{0, 0, 0};
    Vector3D _angularVelocity{0, 0, 0};
};

}  // namespace vox