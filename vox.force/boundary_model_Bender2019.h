//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.force/boundary_model.h"
#include "vox.force/discrete_grid.h"
#include "vox.force/sph_kernels3.h"

namespace vox::flex {
/** \brief The boundary model stores the information required for boundary handling
 * using the approach of Bender et al. 2019 [BKWK19].
 *
 * References:
 * - [BKWK19] Jan Bender, Tassilo Kugelstadt, Marcel Weiler, and Dan Koschier. Volume maps: an implicit boundary
 * representation for SPH. In Proceedings of ACM SIGGRAPH Conference on Motion, Interaction and Games, MIG '19. ACM,
 * 2019. URL: https://dl.acm.org/doi/10.1145/3359566.3360077
 */
class BoundaryModel_Bender2019 : public BoundaryModel {
public:
    BoundaryModel_Bender2019();
    ~BoundaryModel_Bender2019() override;

protected:
    // Density or volume map
    DiscreteGrid *m_map;
    // values required for volume maps
    std::vector<std::vector<double>> m_boundaryVolume;
    std::vector<std::vector<Vector3D>> m_boundaryXj;
    // maxmimal distance of a mesh point to the center of mass (required for CFL)
    double m_maxDist;
    double m_maxVel;

public:
    void initModel(RigidBodyObject *rbo);

    void reset() override;

    DiscreteGrid *getMap() { return m_map; }
    void setMap(DiscreteGrid *map) { m_map = map; }

    [[nodiscard]] double getMaxDist() const { return m_maxDist; }
    void setMaxDist(double val) { m_maxDist = val; }

    [[nodiscard]] double getMaxVel() const { return m_maxVel; }
    void setMaxVel(double val) { m_maxVel = val; }

    [[nodiscard]] inline const double &getBoundaryVolume(const unsigned int fluidIndex, const unsigned int i) const {
        return m_boundaryVolume[fluidIndex][i];
    }

    inline double &getBoundaryVolume(const unsigned int fluidIndex, const unsigned int i) {
        return m_boundaryVolume[fluidIndex][i];
    }

    inline void setBoundaryVolume(const unsigned int fluidIndex, const unsigned int i, const double &val) {
        m_boundaryVolume[fluidIndex][i] = val;
    }

    inline Vector3D &getBoundaryXj(const unsigned int fluidIndex, const unsigned int i) {
        return m_boundaryXj[fluidIndex][i];
    }

    [[nodiscard]] inline const Vector3D &getBoundaryXj(const unsigned int fluidIndex, const unsigned int i) const {
        return m_boundaryXj[fluidIndex][i];
    }

    inline void setBoundaryXj(const unsigned int fluidIndex, const unsigned int i, const Vector3D &val) {
        m_boundaryXj[fluidIndex][i] = val;
    }
};
}  // namespace vox::flex