//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/fluid_model.h"
#include "vox.force/non_pressure_force.h"

namespace vox::flex {
/** \brief Base class for all vorticity methods.
 */
class VorticityBase : public NonPressureForceBase {
protected:
    double m_vorticityCoeff;

    virtual void initParameters();

public:
    static int VORTICITY_COEFFICIENT;

    explicit VorticityBase(FluidModel *model);
    ~VorticityBase() override;
};

}  // namespace vox::flex
