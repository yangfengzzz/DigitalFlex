//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/fluid_model.h"
#include "vox.force/non_pressure_force.h"

namespace vox::flex {
/** \brief Base class for all drag force methods.
 */
class DragBase : public NonPressureForceBase {
protected:
    double m_dragCoefficient;

    virtual void initParameters();

public:
    static int DRAG_COEFFICIENT;

    explicit DragBase(FluidModel *model);
    ~DragBase() override;
};

}  // namespace vox::flex