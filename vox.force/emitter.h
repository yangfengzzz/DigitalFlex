//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.force/fluid_model.h"

namespace vox::flex {
class Emitter {
public:
    Emitter(FluidModel *model,
            const unsigned int width,
            const unsigned int height,
            const Vector3D &pos,
            const Matrix3x3D &rotation,
            const double velocity,
            const unsigned int type = 0);
    virtual ~Emitter();

protected:
    FluidModel *m_model;
    unsigned int m_width;
    unsigned int m_height;
    Vector3D m_x;
    Matrix3x3D m_rotation;
    double m_velocity;
    unsigned int m_type;
    double m_nextEmitTime;
    double m_emitStartTime;
    double m_emitEndTime;
    unsigned int m_emitCounter;
    unsigned int m_objectId;

    static inline bool inBox(const Vector3D &x,
                             const Vector3D &xBox,
                             const Matrix3x3D &rotBox,
                             const Vector3D &scaleBox) {
        const Vector3D xlocal = rotBox.transposed() * (x - xBox);
        // for a box shape, m_scale stores the half-size of the box
        // inside box if closer than half-size on all axes
        return (std::abs(xlocal.x) < scaleBox.x) && (std::abs(xlocal.y) < scaleBox.y) &&
               (std::abs(xlocal.y) < scaleBox.y);
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

public:
    void emitParticles(std::vector<unsigned int> &reusedParticles,
                       unsigned int &indexReuse,
                       unsigned int &numEmittedParticles);
    void emitParticlesCircle(std::vector<unsigned int> &reusedParticles,
                             unsigned int &indexReuse,
                             unsigned int &numEmittedParticles);
    [[nodiscard]] double getNextEmitTime() const { return m_nextEmitTime; }
    void setNextEmitTime(double val) { m_nextEmitTime = val; }
    void setEmitStartTime(double val) {
        m_emitStartTime = val;
        setNextEmitTime(val);
    }
    void setEmitEndTime(double val) { m_emitEndTime = val; }
    static Vector3D getSize(const double width, const double height, const int type);

    void step(std::vector<unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles);
    virtual void reset();

    void saveState(BinaryFileWriter &binWriter);
    void loadState(BinaryFileReader &binReader);

    [[nodiscard]] const Vector3D &getPosition() const { return m_x; }
    void setPosition(const Vector3D &x) { m_x = x; }
    [[nodiscard]] const Matrix3x3D &getRotation() const { return m_rotation; }
    void setRotation(const Matrix3x3D &r) { m_rotation = r; }
    [[nodiscard]] double getVelocity() const { return m_velocity; }
    void setVelocity(const double v) { m_velocity = v; }
    [[nodiscard]] unsigned int getObjectId() const { return m_objectId; }
    void setObjectId(const unsigned int v) { m_objectId = v; }
};
}  // namespace vox::flex