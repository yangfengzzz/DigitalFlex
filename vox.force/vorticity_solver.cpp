//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/vorticity_solver.h"

namespace vox::flex {

int VorticityBase::VORTICITY_COEFFICIENT = -1;

VorticityBase::VorticityBase(FluidModel* model) : NonPressureForceBase(model) {
    m_vorticityCoeff = static_cast<double >(0.01);
}

VorticityBase::~VorticityBase() = default;

void VorticityBase::initParameters() {}

}  // namespace vox::flex