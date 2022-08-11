//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/reflect/parameter_object.h"
#include "vox.force/boundary_model.h"
#include "vox.force/discrete_grid.h"
#include "vox.force/fluid_model.h"

namespace vox::flex {
/** \brief Base class for the simulation methods.
 */
class TimeStep : public utility::ParameterObject {
public:
    static int SOLVER_ITERATIONS;
    static int MIN_ITERATIONS;
    static int MAX_ITERATIONS;
    static int MAX_ERROR;

protected:
    unsigned int m_iterations;
    double m_maxError;
    unsigned int m_minIterations;
    unsigned int m_maxIterations;

    /** Clear accelerations and add gravitation.
     */
    void clearAccelerations(const unsigned int fluidModelIndex);

    virtual void initParameters();

    void approximateNormal(DiscreteGrid *map, const Vector3D &x, Vector3D &n, const unsigned int dim);
    void computeVolumeAndBoundaryX(const unsigned int fluidModelIndex, const unsigned int i, const Vector3D &xi);
    void computeVolumeAndBoundaryX();
    void computeDensityAndGradient(const unsigned int fluidModelIndex, const unsigned int i, const Vector3D &xi);
    void computeDensityAndGradient();

public:
    TimeStep();
    virtual ~TimeStep(void);

    /** Determine densities of all fluid particles.
     */
    void computeDensities(const unsigned int fluidModelIndex);

    virtual void step() = 0;
    virtual void reset();

    virtual void init();
    virtual void resize() = 0;

    virtual void emittedParticles(FluidModel *model, const unsigned int startIndex){};

    virtual void saveState(BinaryFileWriter &binWriter){};
    virtual void loadState(BinaryFileReader &binReader){};

#ifdef USE_PERFORMANCE_OPTIMIZATION
    void precomputeValues();
#endif
};
}  // namespace vox::flex
