// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.force/advection_solver3.h"

using namespace vox;

AdvectionSolver3::AdvectionSolver3() = default;

AdvectionSolver3::~AdvectionSolver3() = default;

void AdvectionSolver3::advect(const CollocatedVectorGrid3& source,
                              const VectorField3& flow,
                              double dt,
                              CollocatedVectorGrid3* target,
                              const ScalarField3& boundarySdf) {}

void AdvectionSolver3::advect(const FaceCenteredGrid3& source,
                              const VectorField3& flow,
                              double dt,
                              FaceCenteredGrid3* target,
                              const ScalarField3& boundarySdf) {}
