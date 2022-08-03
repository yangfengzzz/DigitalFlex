//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "surface_tension_solver.h"

namespace vox ::flex {
int SurfaceTensionBase::SURFACE_TENSION = -1;
int SurfaceTensionBase::SURFACE_TENSION_BOUNDARY = -1;

SurfaceTensionBase::SurfaceTensionBase(FluidModel* model) : NonPressureForceBase(model) {
    m_surfaceTension = 0.05;
    m_surfaceTensionBoundary = 0.01;
}

SurfaceTensionBase::~SurfaceTensionBase() = default;

void SurfaceTensionBase::initParameters() {}

}  // namespace vox::flex
