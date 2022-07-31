// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.force/advection_solver2.h"

using namespace vox;

AdvectionSolver2::AdvectionSolver2() = default;

AdvectionSolver2::~AdvectionSolver2() = default;

void AdvectionSolver2::advect(const CollocatedVectorGrid2& source,
                              const VectorField2& flow,
                              double dt,
                              CollocatedVectorGrid2* target,
                              const ScalarField2& boundarySdf) {}

void AdvectionSolver2::advect(const FaceCenteredGrid2& source,
                              const VectorField2& flow,
                              double dt,
                              FaceCenteredGrid2* target,
                              const ScalarField2& boundarySdf) {}
