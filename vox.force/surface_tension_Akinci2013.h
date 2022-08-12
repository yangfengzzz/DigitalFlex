//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/fluid_model.h"
#include "vox.force/surface_tension_solver.h"

namespace vox::flex {
/** \brief This class implements the surface tension method introduced
 * by Akinci et al. [ATT13].
 *
 * References:
 * - [AAT13] Nadir Akinci, Gizem Akinci, and Matthias Teschner. Versatile surface tension and adhesion for sph fluids.
 * ACM Trans. Graph., 32(6):182:1-182:8, November 2013. URL: http://doi.acm.org/10.1145/2508363.2508395
 */
class SurfaceTension_Akinci2013 : public SurfaceTensionBase {
protected:
    std::vector<Vector3D> m_normals;

public:
    explicit SurfaceTension_Akinci2013(FluidModel *model);
    ~SurfaceTension_Akinci2013() override;

    static NonPressureForceBase *creator(FluidModel *model) { return new SurfaceTension_Akinci2013(model); }

    void step() override;
    void reset() override;

    void computeNormals();

    void performNeighborhoodSearchSort() override;

    inline Vector3D &getNormal(const unsigned int i) { return m_normals[i]; }

    [[nodiscard]] inline const Vector3D &getNormal(const unsigned int i) const { return m_normals[i]; }

    inline void setNormal(const unsigned int i, const Vector3D &val) { m_normals[i] = val; }
};
}  // namespace vox::flex