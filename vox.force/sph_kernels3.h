//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

// Adopted from the sample code of:
// Bart Adams and Martin Wicke,
// "Meshless Approximation Methods and Applications in Physics Based Modeling
// and Animation", Eurographics 2009 Tutorial

#include "vox.base/constants.h"
#include "vox.math/vector3.h"

namespace vox {

//!
//! \brief Standard 3-D SPH kernel function object.
//!
//! \see Müller, Matthias, David Charypar, and Markus Gross.
//!     "Particle-based fluid simulation for interactive applications."
//!     Proceedings of the 2003 ACM SIGGRAPH/Eurographics symposium on Computer
//!     animation. Eurographics Association, 2003.
//!
struct SphStdKernel3 {
    //! Kernel radius.
    double h;

    //! Square of the kernel radius.
    double h2;

    //! Cubic of the kernel radius.
    double h3;

    //! Fifth-power of the kernel radius.
    double h5;

    //! Constructs a kernel object with zero radius.
    SphStdKernel3();

    //! Constructs a kernel object with given radius.
    explicit SphStdKernel3(double kernelRadius);

    //! Copy constructor
    SphStdKernel3(const SphStdKernel3& other);

    //! Returns kernel function value at given distance.
    double operator()(double distance) const;

    //! Returns the first derivative at given distance.
    [[nodiscard]] double firstDerivative(double distance) const;

    //! Returns the gradient at a point.
    [[nodiscard]] Vector3D gradient(const Vector3D& point) const;

    //! Returns the gradient at a point defined by distance and direction.
    [[nodiscard]] Vector3D gradient(double distance, const Vector3D& direction) const;

    //! Returns the second derivative at given distance.
    [[nodiscard]] double secondDerivative(double distance) const;
};

//!
//! \brief Spiky 3-D SPH kernel function object.
//!
//! \see Müller, Matthias, David Charypar, and Markus Gross.
//!     "Particle-based fluid simulation for interactive applications."
//!     Proceedings of the 2003 ACM SIGGRAPH/Eurographics symposium on Computer
//!     animation. Eurographics Association, 2003.
//!
struct SphSpikyKernel3 {
    //! Kernel radius.
    double h;

    //! Square of the kernel radius.
    double h2;

    //! Cubic of the kernel radius.
    double h3;

    //! Fourth-power of the kernel radius.
    double h4;

    //! Fifth-power of the kernel radius.
    double h5;

    //! Constructs a kernel object with zero radius.
    SphSpikyKernel3();

    //! Constructs a kernel object with given radius.
    explicit SphSpikyKernel3(double h_);

    //! Copy constructor
    SphSpikyKernel3(const SphSpikyKernel3& other);

    //! Returns kernel function value at given distance.
    double operator()(double distance) const;

    //! Returns the first derivative at given distance.
    [[nodiscard]] double firstDerivative(double distance) const;

    //! Returns the gradient at a point.
    [[nodiscard]] Vector3D gradient(const Vector3D& point) const;

    //! Returns the gradient at a point defined by distance and direction.
    [[nodiscard]] Vector3D gradient(double distance, const Vector3D& direction) const;

    //! Returns the second derivative at given distance.
    [[nodiscard]] double secondDerivative(double distance) const;
};

}  // namespace vox

#include "vox.force/sph_kernels3-inl.h"