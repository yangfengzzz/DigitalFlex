//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/fluid_model.h"
#include "vox.force/non_pressure_force.h"

namespace vox::flex {
/** \brief Base class for all surface tension methods.
 */
class SurfaceTensionBase : public NonPressureForceBase {
protected:
    double m_surfaceTension;
    double m_surfaceTensionBoundary;

    void initParameters() override;

public:
    static int SURFACE_TENSION;
    static int SURFACE_TENSION_BOUNDARY;

    explicit SurfaceTensionBase(FluidModel *model);
    ~SurfaceTensionBase() override;
};

}  // namespace vox::flex