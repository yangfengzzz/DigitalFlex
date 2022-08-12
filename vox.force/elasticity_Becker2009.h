//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/elasticity_solver.h"
#include "vox.force/fluid_model.h"

namespace vox::flex {
/** \brief This class implements the corotated SPH method for deformable solids introduced
 * by Becker et al. [BIT09].
 *
 * References:
 * - [BIT09] Markus Becker, Markus Ihmsen, and Matthias Teschner. Corotated SPH for deformable solids. In Proceedings of
 * Eurographics Conference on Natural Phenomena, 27-34. 2009. URL:
 * http://dx.doi.org/10.2312EG/DL/conf/EG2009/nph/027-034
 */
class Elasticity_Becker2009 : public ElasticityBase {
protected:
    // initial particle indices, used to access their original positions
    std::vector<unsigned int> m_current_to_initial_index;
    std::vector<unsigned int> m_initial_to_current_index;
    // initial particle neighborhood
    std::vector<std::vector<unsigned int>> m_initialNeighbors;
    // volumes in rest configuration
    std::vector<double> m_restVolumes;
    std::vector<Matrix3x3D> m_rotations;
    std::vector<Vector<double, 6>> m_stress;
    std::vector<Matrix3x3D> m_F;
    double m_alpha;

    void initValues();
    void computeRotations();
    void computeStress();
    void computeForces();

    void initParameters() override;

    //////////////////////////////////////////////////////////////////////////
    // multiplication of symmetric matrix, represented by a 6D vector, and a
    // 3D vector
    //////////////////////////////////////////////////////////////////////////
    static inline void symMatTimesVec(const Vector<double, 6> &M, const Vector3D &v, Vector3D &res) {
        res[0] = M[0] * v[0] + M[3] * v[1] + M[4] * v[2];
        res[1] = M[3] * v[0] + M[1] * v[1] + M[5] * v[2];
        res[2] = M[4] * v[0] + M[5] * v[1] + M[2] * v[2];
    }

public:
    static int ALPHA;

    explicit Elasticity_Becker2009(FluidModel *model);
    ~Elasticity_Becker2009() override;

    static NonPressureForceBase *creator(FluidModel *model) { return new Elasticity_Becker2009(model); }

    void step() override;
    void reset() override;
    void performNeighborhoodSearchSort() override;

    void saveState(BinaryFileWriter &binWriter) override;
    void loadState(BinaryFileReader &binReader) override;
};
}  // namespace vox::flex
