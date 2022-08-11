//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.force/boundary_model.h"
#include "vox.force/sph_kernels3.h"

namespace vox::flex {
/** \brief The boundary model stores the information required for boundary handling
 * using the approach of Akinci et al. 2012 [AIA+12].
 *
 * References:
 * - [AIA+12] Nadir Akinci, Markus Ihmsen, Gizem Akinci, Barbara Solenthaler, and Matthias Teschner. Versatile
 * rigid-fluid coupling for incompressible SPH. ACM Trans. Graph., 31(4):62:1-62:8, July 2012. URL:
 * http://doi.acm.org/10.1145/2185520.2185558
 */
class BoundaryModel_Akinci2012 : public BoundaryModel {
public:
    BoundaryModel_Akinci2012();
    ~BoundaryModel_Akinci2012() override;

protected:
    bool m_sorted;
    unsigned int m_pointSetIndex;

    // values required for Akinci 2012 boundary handling
    std::vector<Vector3D> m_x0;
    std::vector<Vector3D> m_x;
    std::vector<Vector3D> m_v;
    std::vector<double> m_V;

public:
    [[nodiscard]] unsigned int numberOfParticles() const { return static_cast<unsigned int>(m_x.size()); }
    [[nodiscard]] unsigned int getPointSetIndex() const { return m_pointSetIndex; }
    [[nodiscard]] bool isSorted() const { return m_sorted; }

    void computeBoundaryVolume();
    void resize(unsigned int numBoundaryParticles);

    void reset() override;

    void performNeighborhoodSearchSort() override;

    void saveState(BinaryFileWriter &binWriter) override;
    void loadState(BinaryFileReader &binReader) override;

    void initModel(RigidBodyObject *rbo, unsigned int numBoundaryParticles, Vector3D *boundaryParticles);

    inline Vector3D &getPosition0(const unsigned int i) { return m_x0[i]; }

    [[nodiscard]] inline const Vector3D &getPosition0(const unsigned int i) const { return m_x0[i]; }

    inline void setPosition0(const unsigned int i, const Vector3D &pos) { m_x0[i] = pos; }

    inline Vector3D &getPosition(const unsigned int i) { return m_x[i]; }

    [[nodiscard]] inline const Vector3D &getPosition(const unsigned int i) const { return m_x[i]; }

    inline void setPosition(const unsigned int i, const Vector3D &pos) { m_x[i] = pos; }

    inline Vector3D &getVelocity(const unsigned int i) { return m_v[i]; }

    [[nodiscard]] inline const Vector3D &getVelocity(const unsigned int i) const { return m_v[i]; }

    inline void setVelocity(const unsigned int i, const Vector3D &vel) { m_v[i] = vel; }

    [[nodiscard]] inline const double &getVolume(const unsigned int i) const { return m_V[i]; }

    inline double &getVolume(const unsigned int i) { return m_V[i]; }

    inline void setVolume(const unsigned int i, const double &val) { m_V[i] = val; }
};
}  // namespace vox::flex