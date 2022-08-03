//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.force/animation_field.h"

namespace vox::flex {
class FluidModel;

class AnimationFieldSystem {
public:
    AnimationFieldSystem();
    virtual ~AnimationFieldSystem();

protected:
    std::vector<AnimationField *> m_fields;

public:
    void addAnimationField(const std::string &particleFieldName,
                           const Vector3D &pos,
                           const Matrix3x3D &rotation,
                           const Vector3D &scale,
                           const std::string expression[3],
                           unsigned int type);
    [[nodiscard]] unsigned int numAnimationFields() const { return static_cast<unsigned int>(m_fields.size()); }
    std::vector<AnimationField *> &getAnimationFields() { return m_fields; }

    void step();
    void reset();
};
}  // namespace vox::flex