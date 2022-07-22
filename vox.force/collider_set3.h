//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.force/collider3.h"
#include "vox.geometry/surface_set3.h"

namespace vox {

//! Collection of 3-D colliders
class ColliderSet3 final : public Collider3 {
public:
    class Builder;

    //! Default constructor.
    ColliderSet3();

    //! Constructs with other colliders.
    explicit ColliderSet3(const std::vector<Collider3Ptr>& others);

    //! Returns the velocity of the collider at given \p point.
    [[nodiscard]] Vector3D velocityAt(const Point3D& point) const override;

    //! Adds a collider to the set.
    void addCollider(const Collider3Ptr& collider);

    //! Returns number of colliders.
    [[nodiscard]] size_t numberOfColliders() const;

    //! Returns collider at index \p i.
    [[nodiscard]] Collider3Ptr collider(size_t i) const;

    //! Returns builder fox ColliderSet3.
    static Builder builder();

private:
    std::vector<Collider3Ptr> _colliders;
};

//! Shared pointer for the ColliderSet3 type.
typedef std::shared_ptr<ColliderSet3> ColliderSet3Ptr;

//!
//! \brief Front-end to create ColliderSet3 objects step by step.
//!
class ColliderSet3::Builder final {
public:
    //! Returns builder with other colliders.
    Builder& withColliders(const std::vector<Collider3Ptr>& others);

    //! Builds ColliderSet3.
    [[nodiscard]] ColliderSet3 build() const;

    //! Builds shared pointer of ColliderSet3 instance.
    [[nodiscard]] ColliderSet3Ptr makeShared() const;

private:
    std::vector<Collider3Ptr> _colliders;
};

}  // namespace vox