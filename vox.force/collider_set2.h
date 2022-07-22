//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.force/collider2.h"
#include "vox.geometry/surface_set2.h"

namespace vox {

//! Collection of 2-D colliders
class ColliderSet2 final : public Collider2 {
public:
    class Builder;

    //! Default constructor.
    ColliderSet2();

    //! Constructs with other colliders.
    explicit ColliderSet2(const std::vector<Collider2Ptr>& others);

    //! Returns the velocity of the collider at given \p point.
    [[nodiscard]] Vector2D velocityAt(const Point2D& point) const override;

    //! Adds a collider to the set.
    void addCollider(const Collider2Ptr& collider);

    //! Returns number of colliders.
    [[nodiscard]] size_t numberOfColliders() const;

    //! Returns collider at index \p i.
    [[nodiscard]] Collider2Ptr collider(size_t i) const;

    //! Returns builder fox ColliderSet2.
    static Builder builder();

private:
    std::vector<Collider2Ptr> _colliders;
};

//! Shared pointer for the ColliderSet2 type.
typedef std::shared_ptr<ColliderSet2> ColliderSet2Ptr;

//!
//! \brief Front-end to create ColliderSet2 objects step by step.
//!
class ColliderSet2::Builder final {
public:
    //! Returns builder with other colliders.
    Builder& withColliders(const std::vector<Collider2Ptr>& others);

    //! Builds ColliderSet2.
    [[nodiscard]] ColliderSet2 build() const;

    //! Builds shared pointer of ColliderSet2 instance.
    [[nodiscard]] ColliderSet2Ptr makeShared() const;

private:
    std::vector<Collider2Ptr> _colliders;
};

}  // namespace vox