//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/animation_field_system.h"
#include "vox.force/fluid_model.h"

namespace vox::flex {

AnimationFieldSystem::AnimationFieldSystem() : m_fields() {}

AnimationFieldSystem::~AnimationFieldSystem() {
    for (auto & m_field : m_fields) delete m_field;
}

void AnimationFieldSystem::step() {
    for (auto & m_field : m_fields) {
        m_field->step();
    }
}

void AnimationFieldSystem::reset() {
    for (auto & m_field : m_fields) {
        m_field->reset();
    }
}

void AnimationFieldSystem::addAnimationField(const std::string &particleFieldName,
                                             const Vector3D &pos,
                                             const Matrix3x3D &rotation,
                                             const Vector3D &scale,
                                             const std::string expression[3],
                                             const unsigned int type) {
    m_fields.push_back(new AnimationField(particleFieldName, pos, rotation, scale, expression, type));
}
}  // namespace vox::flex