//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/fluid_model.h"
#include "vox.force/non_pressure_force.h"

namespace vox::flex {
/**
 * \brief Base class for all viscosity methods.
 */
class ViscosityBase : public NonPressureForceBase {
protected:
    double m_viscosity;

    virtual void initParameters();

public:
    static int VISCOSITY_COEFFICIENT;

    explicit ViscosityBase(FluidModel *model);
    ~ViscosityBase() override;
};

}  // namespace vox::flex