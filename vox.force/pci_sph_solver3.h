//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/sph_solver3.h"

namespace vox {

//!
//! \brief 3-D PCISPH solver.
//!
//! This class implements 3-D predictive-corrective SPH solver. The main
//! pressure solver is based on Solenthaler and Pajarola's 2009 SIGGRAPH paper.
//!
//! \see Solenthaler and Pajarola, Predictive-corrective incompressible SPH,
//!      ACM transactions on graphics (TOG). Vol. 28. No. 3. ACM, 2009.
//!
class PciSphSolver3 : public SphSolver3 {
public:
    class Builder;

    //! Constructs a solver with empty particle set.
    PciSphSolver3();

    //! Constructs a solver with target density, spacing, and relative kernel
    //! radius.
    PciSphSolver3(double targetDensity, double targetSpacing, double relativeKernelRadius);

    ~PciSphSolver3() override;

    //! Returns max allowed density error ratio.
    [[nodiscard]] double maxDensityErrorRatio() const;

    //!
    //! \brief Sets max allowed density error ratio.
    //!
    //! This function sets the max allowed density error ratio during the PCISPH
    //! iteration. Default is 0.01 (1%). The input value should be positive.
    //!
    void setMaxDensityErrorRatio(double ratio);

    //! Returns max number of iterations.
    [[nodiscard]] unsigned int maxNumberOfIterations() const;

    //!
    //! \brief Sets max number of PCISPH iterations.
    //!
    //! This function sets the max number of PCISPH iterations. Default is 5.
    //!
    void setMaxNumberOfIterations(unsigned int n);

    //! Returns builder fox PciSphSolver3.
    static Builder builder();

protected:
    //! Accumulates the pressure force to the forces array in the particle
    //! system.
    void accumulatePressureForce(double timeIntervalInSeconds) override;

    //! Performs pre-processing step before the simulation.
    void onBeginAdvanceTimeStep(double timeStepInSeconds) override;

private:
    double _maxDensityErrorRatio = 0.01;
    unsigned int _maxNumberOfIterations = 5;

    Array1<Point3D> _tempPositions;
    ParticleSystemData3::VectorData _tempVelocities;
    ParticleSystemData3::VectorData _pressureForces;
    ParticleSystemData3::ScalarData _densityErrors;

    double computeDelta(double timeStepInSeconds);
    double computeBeta(double timeStepInSeconds);
};

//! Shared pointer type for the PciSphSolver3.
typedef std::shared_ptr<PciSphSolver3> PciSphSolver3Ptr;

//!
//! \brief Front-end to create PciSphSolver3 objects step by step.
//!
class PciSphSolver3::Builder final : public SphSolverBuilderBase3<PciSphSolver3::Builder> {
public:
    //! Builds PciSphSolver3.
    [[nodiscard]] PciSphSolver3 build() const;

    //! Builds shared pointer of PciSphSolver3 instance.
    [[nodiscard]] PciSphSolver3Ptr makeShared() const;
};

}  // namespace vox