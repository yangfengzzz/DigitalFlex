//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/fluid_model.h"
#include "vox.force/vorticity_solver.h"

namespace vox::flex {
/** \brief This class implements the vorticity confinement method introduced
 * by Macklin and Mueller [MM13].
 *
 * References:
 * - [MM13] Miles Macklin and Matthias Müller. Position based fluids. ACM Trans. Graph., 32(4):104:1-104:12, July 2013.
 * URL: http://doi.acm.org/10.1145/2461912.2461984
 */
class VorticityConfinement : public VorticityBase {
protected:
    std::vector<Vector3D> m_omega;
    std::vector<double> m_normOmega;

public:
    explicit VorticityConfinement(FluidModel* model);
    ~VorticityConfinement() override;

    static NonPressureForceBase* creator(FluidModel* model) { return new VorticityConfinement(model); }

    void step() override;
    void reset() override;

    void performNeighborhoodSearchSort() override;

    [[nodiscard]] inline const Vector3D& getAngularVelocity(const unsigned int i) const { return m_omega[i]; }

    inline Vector3D& getAngularVelocity(const unsigned int i) { return m_omega[i]; }

    inline void setAngularVelocity(const unsigned int i, const Vector3D& val) { m_omega[i] = val; }
};
}  // namespace vox::flex
