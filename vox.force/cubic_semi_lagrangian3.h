//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/semi_lagrangian3.h"

namespace vox {

//!
//! \brief Implementation of 3-D cubic semi-Lagrangian advection solver.
//!
//! This class implements 3rd-order cubic 3-D semi-Lagrangian advection solver.
//!
class CubicSemiLagrangian3 final : public SemiLagrangian3 {
public:
    CubicSemiLagrangian3();

protected:
    //!
    //! \brief Returns spatial interpolation function object for given scalar
    //! grid.
    //!
    //! This function overrides the original function with cubic interpolation.
    //!
    [[nodiscard]] std::function<double(const Point3D&)> getScalarSamplerFunc(const ScalarGrid3& source) const override;

    //!
    //! \brief Returns spatial interpolation function object for given
    //! collocated vector grid.
    //!
    //! This function overrides the original function with cubic interpolation.
    //!
    [[nodiscard]] std::function<Vector3D(const Point3D&)> getVectorSamplerFunc(
            const CollocatedVectorGrid3& source) const override;

    //!
    //! \brief Returns spatial interpolation function object for given
    //! face-centered vector grid.
    //!
    //! This function overrides the original function with cubic interpolation.
    //!
    [[nodiscard]] std::function<Vector3D(const Point3D&)> getVectorSamplerFunc(
            const FaceCenteredGrid3& source) const override;
};

typedef std::shared_ptr<CubicSemiLagrangian3> CubicSemiLagrangian3Ptr;

}  // namespace vox
