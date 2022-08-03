//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/fluid_model.h"
#include "vox.force/non_pressure_force.h"

namespace vox::flex {
/** \brief Base class for all elasticity methods.
 */
class ElasticityBase : public NonPressureForceBase {
protected:
    double m_youngsModulus;
    double m_poissonRatio;
    Vector3D m_fixedBoxMin;
    Vector3D m_fixedBoxMax;

    virtual void initParameters();
    void determineFixedParticles();

public:
    static int YOUNGS_MODULUS;
    static int POISSON_RATIO;
    static int FIXED_BOX_MIN;
    static int FIXED_BOX_MAX;

    explicit ElasticityBase(FluidModel *model);
    ~ElasticityBase() override;
};
}  // namespace vox::flex