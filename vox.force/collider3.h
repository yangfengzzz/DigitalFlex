//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <functional>

#include "vox.geometry/surface3.h"

namespace vox {

//!
//! \brief Abstract base class for generic collider object.
//!
//! This class contains basic interfaces for colliders. Most of the
//! functionalities are implemented within this class, except the member
//! function Collider3::velocityAt. This class also let the subclasses to
//! provide a Surface3 instance to define collider surface using
//! Collider3::setSurface function.
//!
class Collider3 {
public:
    //!
    //! \brief Callback function type for update calls.
    //!
    //! This type of callback function will take the collider pointer, current
    //! time, and time interval in seconds.
    //!
    typedef std::function<void(Collider3*, double, double)> OnBeginUpdateCallback;

    //! Default constructor.
    Collider3();

    //! Default destructor.
    virtual ~Collider3();

    //! Returns the velocity of the collider at given \p point.
    [[nodiscard]] virtual Vector3D velocityAt(const Point3D& point) const = 0;

    //!
    //! Resolves collision for given point.
    //!
    //! \param radius Radius of the colliding point.
    //! \param restitutionCoefficient Defines the restitution effect.
    //! \param position Input and output position of the point.
    //! \param velocity Input and output velocity of the point.
    //!
    void resolveCollision(double radius, double restitutionCoefficient, Point3D* position, Vector3D* velocity);

    //! Returns friction coefficent.
    [[nodiscard]] double frictionCoefficient() const;

    //!
    //! \brief Sets the friction coefficient.
    //!
    //! This function assigns the friction coefficient to the collider. Any
    //! negative inputs will be clamped to zero.
    //!
    void setFrictionCoefficient(double newFrictionCoeffient);

    //! Returns the surface instance.
    [[nodiscard]] const Surface3Ptr& surface() const;

    //! Updates the collider state.
    void update(double currentTimeInSeconds, double timeIntervalInSeconds);

    //!
    //! \brief      Sets the callback function to be called when
    //!             Collider2::update function is invoked.
    //!
    //! The callback function takes current simulation time in seconds unit. Use
    //! this callback to track any motion or state changes related to this
    //! collider.
    //!
    //! \param[in]  callback The callback function.
    //!
    void setOnBeginUpdateCallback(const OnBeginUpdateCallback& callback);

protected:
    //! Internal query result structure.
    struct ColliderQueryResult final {
        double distance;
        Point3D point;
        Vector3D normal;
        Vector3D velocity;
    };

    //! Assigns the surface instance from the subclass.
    void setSurface(const Surface3Ptr& newSurface);

    //! Outputs closest point's information.
    void getClosestPoint(const Surface3Ptr& surface, const Point3D& queryPoint, ColliderQueryResult* result) const;

    //! Returns true if given point is in the opposite side of the surface.
    bool isPenetrating(const ColliderQueryResult& colliderPoint, const Point3D& position, double radius);

private:
    Surface3Ptr _surface;
    double _frictionCoeffient = 0.0;
    OnBeginUpdateCallback _onUpdateCallback;
};

//! Shared pointer type for the Collider2.
typedef std::shared_ptr<Collider3> Collider3Ptr;

}  // namespace vox