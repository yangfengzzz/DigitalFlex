//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/viscosity_solver.h"

namespace vox::flex {

int ViscosityBase::VISCOSITY_COEFFICIENT = -1;

ViscosityBase::ViscosityBase(FluidModel* model) : NonPressureForceBase(model) { m_viscosity = 0.01; }

ViscosityBase::~ViscosityBase() = default;

void ViscosityBase::initParameters() {}

}  // namespace vox::flex
