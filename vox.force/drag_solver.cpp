//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/drag_solver.h"

namespace vox::flex {

int DragBase::DRAG_COEFFICIENT = -1;

DragBase::DragBase(FluidModel* model) : NonPressureForceBase(model) { m_dragCoefficient = static_cast<double>(0.01); }

DragBase::~DragBase() = default;

void DragBase::initParameters() {}

}  // namespace vox::flex
