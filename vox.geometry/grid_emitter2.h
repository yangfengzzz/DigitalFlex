//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <utility>
#include <vector>

#include "vox.geometry/implicit_surface2.h"
#include "vox.geometry/scalar_grid2.h"

namespace vox {

//!
//! \brief Abstract base class for 2-D grid-based emitters.
//!
class GridEmitter2 {
public:
    //!
    //! \brief Callback function type for update calls.
    //!
    //! This type of callback function will take the current time and time
    //! interval in seconds.
    //!
    typedef std::function<void(double, double)> OnBeginUpdateCallback;

    //! Constructs an emitter.
    GridEmitter2();

    //! Destructor.
    virtual ~GridEmitter2();

    //! Updates the emitter state from \p currentTimeInSeconds to the following
    //! time-step.
    void update(double currentTimeInSeconds, double timeIntervalInSeconds);

    //! Returns true if the emitter is enabled.
    [[nodiscard]] bool isEnabled() const;

    //! Sets true/false to enable/disable the emitter.
    void setIsEnabled(bool enabled);

    //!
    //! \brief      Sets the callback function to be called when
    //!             GridEmitter2::update function is invoked.
    //!
    //! The callback function takes current simulation time in seconds unit. Use
    //! this callback to track any motion or state changes related to this
    //! emitter.
    //!
    //! \param[in]  callback The callback function.
    //!
    void setOnBeginUpdateCallback(const OnBeginUpdateCallback& callback);

protected:
    virtual void onUpdate(double currentTimeInSeconds, double timeIntervalInSeconds) = 0;

private:
    bool _isEnabled = true;
    OnBeginUpdateCallback _onBeginUpdateCallback;
};

//! Shared pointer type for the GridEmitter2.
typedef std::shared_ptr<GridEmitter2> GridEmitter2Ptr;

}  // namespace vox
