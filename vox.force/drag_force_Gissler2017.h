//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/drag_solver.h"
#include "vox.force/fluid_model.h"

namespace vox::flex {
/** \brief This class implements the drag force computation introduced
 * by Gissler et al. [GBP+17].
 *
 * References:
 * - [GPB+17] Christoph Gissler, Stefan Band, Andreas Peer, Markus Ihmsen, and Matthias Teschner. Approximate air-fluid
 * interactions for SPH. In Virtual Reality Interactions and Physical Simulations, 1-10. April 2017. URL:
 * http://dx.doi.org/10.2312/vriphys.20171081
 */
class DragForce_Gissler2017 : public DragBase {
protected:
    const double rho_a = static_cast<double>(1.2041);
    const double sigma = static_cast<double>(0.0724);
    const double mu_l = static_cast<double>(0.00102);
    const double C_F = static_cast<double>(1.0 / 3.0);
    const double C_k = static_cast<double>(8.0);
    const double C_d = static_cast<double>(5.0);
    const double C_b = static_cast<double>(0.5);
    const double mu_a = static_cast<double>(0.00001845);

public:
    explicit DragForce_Gissler2017(FluidModel* model);
    ~DragForce_Gissler2017() override;

    static NonPressureForceBase* creator(FluidModel* model) { return new DragForce_Gissler2017(model); }

    void step() override;
    void reset() override;
};

}  // namespace vox::flex