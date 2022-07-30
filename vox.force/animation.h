//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <memory>

#include "vox.base/macros.h"
#include "vox.geometry/frame.h"

namespace vox {
//!
//! \brief Abstract base class for animation-related class.
//!
//! This class represents the animation logic in very abstract level.
//! Generally animation is a function of time and/or its previous state.
//! This base class provides a virtual function update() which can be
//! override by its sub-classes to implement their own state update logic.
//!
class Animation {
public:
    Animation();

    virtual ~Animation();

    //!
    //! \brief Updates animation state for given \p frame.
    //!
    //! This function updates animation state by calling Animation::onUpdate
    //! function.
    //!
    void update(const Frame& frame);

protected:
    //!
    //! \brief The implementation of this function should update the animation
    //!     state for given Frame instance \p frame.
    //!
    //! This function is called from Animation::update when state of this class
    //! instance needs to be updated. Thus, the inherited class should override
    //! this function and implement its logic for updating the animation state.
    //!
    virtual void onUpdate(const Frame& frame) = 0;
};

//! Shared pointer for the Animation type.
typedef std::shared_ptr<Animation> AnimationPtr;

}  // namespace vox