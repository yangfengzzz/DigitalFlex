//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.force/fluid_model.h"

namespace vox::flex {
class AnimationField {
public:
    AnimationField(const std::string &particleFieldName,
                   const Vector3D &pos,
                   const Matrix3x3D &rotation,
                   const Vector3D &scale,
                   const std::string expression[3],
                   const unsigned int type = 0);
    virtual ~AnimationField();

protected:
    std::string m_particleFieldName;
    Vector3D m_x;
    Matrix3x3D m_rotation;
    Vector3D m_scale;
    std::string m_expression[3];
    unsigned int m_type;
    double m_startTime;
    double m_endTime;

    static inline bool inBox(const Vector3D &x,
                             const Vector3D &xBox,
                             const Matrix3x3D &rotBox,
                             const Vector3D &scaleBox) {
        const Vector3D xlocal = rotBox.transposed() * (x - xBox);
        // for a box shape, m_scale stores the half-size of the box
        // inside box if closer than half-size on all axes
        return (std::abs(xlocal.x) < scaleBox.x) && (std::abs(xlocal.y) < scaleBox.y) &&
               (std::abs(xlocal.z) < scaleBox.z);
    }

    static inline bool inCylinder(
            const Vector3D &x, const Vector3D &xCyl, const Matrix3x3D &rotCyl, const double h, const double r2) {
        const Vector3D xlocal = rotCyl.transposed() * (x - xCyl);
        // inside cylinder if distance to x-axis is less than r
        // and projection on x-axis is between 0 and h
        const double proj = xlocal.x;
        const double d2 = Vector2D(xlocal.y, xlocal.z).lengthSquared();
        const double hHalf = static_cast<double>(0.5) * h;
        return (proj > -hHalf) && (proj < hHalf) && (d2 < r2);
    }

    static inline bool inSphere(const Vector3D &x, const Vector3D &pos, const Matrix3x3D &rot, const double radius) {
        const Vector3D xlocal = rot.transposed() * (x - pos);
        return xlocal.lengthSquared() < radius * radius;
    }

    static inline bool inShape(
            const int type, const Vector3D &x, const Vector3D &pos, const Matrix3x3D &rot, const Vector3D &scale) {
        if (type == 1)
            return inSphere(x, pos, rot, scale[0]);
        else if (type == 2) {
            const double h = scale[0];
            const double r = scale[1];
            return inCylinder(x, pos, rot, h, r * r);
        } else
            return inBox(x, pos, rot, 0.5 * scale);
    }

public:
    void setStartTime(double val) { m_startTime = val; }
    void setEndTime(double val) { m_endTime = val; }

    void step();
    virtual void reset();
};

}  // namespace vox::flex