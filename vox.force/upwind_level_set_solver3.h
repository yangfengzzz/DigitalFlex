//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/iterative_level_set_solver3.h"

namespace vox {

//! Three-dimensional first-order upwind-based iterative level set solver.
class UpwindLevelSetSolver3 final : public IterativeLevelSetSolver3 {
public:
    //! Default constructor.
    UpwindLevelSetSolver3();

protected:
    //! Computes the derivatives for given grid point.
    void getDerivatives(ConstArrayAccessor3<double> grid,
                        const Vector3D& gridSpacing,
                        size_t i,
                        size_t j,
                        size_t k,
                        std::array<double, 2>* dx,
                        std::array<double, 2>* dy,
                        std::array<double, 2>* dz) const override;
};

typedef std::shared_ptr<UpwindLevelSetSolver3> UpwindLevelSetSolver3Ptr;

}  // namespace vox