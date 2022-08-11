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
 * using the approach of Koschier and Bender 2017 [KB17].
 *
 * References:
 * - [KB17] Dan Koschier and Jan Bender. Density maps for improved SPH boundary handling. In ACM SIGGRAPH/Eurographics
 * Symposium on Computer Animation, 1-10. July 2017. URL: http://dx.doi.org/10.1145/3099564.3099565
 */
class BoundaryModel_Koschier2017 : public BoundaryModel {
public:
    BoundaryModel_Koschier2017();
    ~BoundaryModel_Koschier2017() override;

protected:
    // Density map
    DiscreteGrid *m_map;
    // values required for density maps
    std::vector<std::vector<double>> m_boundaryDensity;
    std::vector<std::vector<Vector3D>> m_boundaryDensityGradient;
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

    [[nodiscard]] inline const double &getBoundaryDensity(const unsigned int fluidIndex, const unsigned int i) const {
        return m_boundaryDensity[fluidIndex][i];
    }

    inline double &getBoundaryDensity(const unsigned int fluidIndex, const unsigned int i) {
        return m_boundaryDensity[fluidIndex][i];
    }

    inline void setBoundaryDensity(const unsigned int fluidIndex, const unsigned int i, const double &val) {
        m_boundaryDensity[fluidIndex][i] = val;
    }

    inline Vector3D &getBoundaryDensityGradient(const unsigned int fluidIndex, const unsigned int i) {
        return m_boundaryDensityGradient[fluidIndex][i];
    }

    [[nodiscard]] inline const Vector3D &getBoundaryDensityGradient(const unsigned int fluidIndex,
                                                                    const unsigned int i) const {
        return m_boundaryDensityGradient[fluidIndex][i];
    }

    inline void setBoundaryDensityGradient(const unsigned int fluidIndex, const unsigned int i, const Vector3D &val) {
        m_boundaryDensityGradient[fluidIndex][i] = val;
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