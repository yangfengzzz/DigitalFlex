//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.math/matrix3x3.h"
#include "vox.math/quaternion.h"
#include "vox.math/vector3.h"

namespace vox::flex {
/** \brief Base class for rigid body objects.
 */
class RigidBodyObject {
protected:
    bool m_isAnimated;

public:
    RigidBodyObject() { m_isAnimated = false; }
    virtual ~RigidBodyObject() = default;

    [[nodiscard]] virtual bool isDynamic() const = 0;

    [[nodiscard]] bool isAnimated() const { return m_isAnimated; }

    virtual void setIsAnimated(const bool b) { m_isAnimated = true; }

    [[nodiscard]] virtual double const getMass() const = 0;

    [[nodiscard]] virtual Vector3D const &getPosition() const = 0;

    virtual void setPosition(const Vector3D &x) = 0;

    [[nodiscard]] virtual Vector3D getWorldSpacePosition() const = 0;

    [[nodiscard]] virtual Vector3D const &getVelocity() const = 0;

    virtual void setVelocity(const Vector3D &v) = 0;

    [[nodiscard]] virtual QuaternionD const &getRotation() const = 0;

    virtual void setRotation(const QuaternionD &q) = 0;

    [[nodiscard]] virtual Matrix3x3D getWorldSpaceRotation() const = 0;

    [[nodiscard]] virtual Vector3D const &getAngularVelocity() const = 0;

    virtual void setAngularVelocity(const Vector3D &v) = 0;

    virtual void addForce(const Vector3D &f) = 0;

    virtual void addTorque(const Vector3D &t) = 0;

    virtual void updateMeshTransformation() = 0;

    [[nodiscard]] virtual const std::vector<Vector3D> &getVertices() const = 0;

    [[nodiscard]] virtual const std::vector<Vector3D> &getVertexNormals() const = 0;

    [[nodiscard]] virtual const std::vector<unsigned int> &getFaces() const = 0;
};
}  // namespace vox::flex