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
    SphStdKernel3(const SphStdKernel3 &other);

    //! Returns kernel function value at given distance.
    double operator()(double distance) const;

    //! Returns the first derivative at given distance.
    [[nodiscard]] double firstDerivative(double distance) const;

    //! Returns the gradient at a point.
    [[nodiscard]] Vector3D gradient(const Vector3D &point) const;

    //! Returns the gradient at a point defined by distance and direction.
    [[nodiscard]] Vector3D gradient(double distance, const Vector3D &direction) const;

    //! Returns the second derivative at given distance.
    [[nodiscard]] double secondDerivative(double distance) const;
};

// MARK: - SphSpikyKernel3
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
    SphSpikyKernel3(const SphSpikyKernel3 &other);

    //! Returns kernel function value at given distance.
    double operator()(double distance) const;

    //! Returns the first derivative at given distance.
    [[nodiscard]] double firstDerivative(double distance) const;

    //! Returns the gradient at a point.
    [[nodiscard]] Vector3D gradient(const Vector3D &point) const;

    //! Returns the gradient at a point defined by distance and direction.
    [[nodiscard]] Vector3D gradient(double distance, const Vector3D &direction) const;

    //! Returns the second derivative at given distance.
    [[nodiscard]] double secondDerivative(double distance) const;
};

// MARK: - CubicKernel
/**
 * \brief Cubic spline kernel.
 */
class CubicKernel {
protected:
    static double m_radius;
    static double m_k;
    static double m_l;
    static double m_W_zero;

public:
    static double getRadius() { return m_radius; }
    static void setRadius(double val) {
        m_radius = val;
        const auto pi = static_cast<double>(M_PI);

        const double h3 = m_radius * m_radius * m_radius;
        m_k = static_cast<double>(8.0) / (pi * h3);
        m_l = static_cast<double>(48.0) / (pi * h3);
        m_W_zero = W(Vector3D());
    }

public:
    static double W(const double r) {
        double res = 0.0;
        const double q = r / m_radius;
        if (q <= 1.0) {
            if (q <= 0.5) {
                const double q2 = q * q;
                const double q3 = q2 * q;
                res = m_k * (static_cast<double>(6.0) * q3 - static_cast<double>(6.0) * q2 + static_cast<double>(1.0));
            } else {
                res = m_k * (static_cast<double>(2.0) * pow(static_cast<double>(1.0) - q, static_cast<double>(3.0)));
            }
        }
        return res;
    }

    static double W(const Vector3D &r) { return W(r.length()); }

    static Vector3D gradW(const Vector3D &r) {
        Vector3D res;
        const double rl = r.length();
        const double q = rl / m_radius;
        if ((rl > 1.0e-5) && (q <= 1.0)) {
            const Vector3D gradq = r * (static_cast<double>(1.0) / (rl * m_radius));
            if (q <= 0.5) {
                res = m_l * q * ((double)3.0 * q - static_cast<double>(2.0)) * gradq;
            } else {
                const double factor = static_cast<double>(1.0) - q;
                res = m_l * (-factor * factor) * gradq;
            }
        } else
            res.setZero();

        return res;
    }

    static double W_zero() { return m_W_zero; }
};

// MARK: - Poly6Kernel
/**
 * \brief Poly6 kernel.
 */
class Poly6Kernel {
protected:
    static double m_radius;
    static double m_k;
    static double m_l;
    static double m_m;
    static double m_W_zero;

public:
    static double getRadius() { return m_radius; }
    static void setRadius(double val) {
        m_radius = val;
        const auto pi = static_cast<double>(M_PI);
        m_k = static_cast<double>(315.0) / (static_cast<double>(64.0) * pi * pow(m_radius, static_cast<double>(9.0)));
        m_l = -static_cast<double>(945.0) / (static_cast<double>(32.0) * pi * pow(m_radius, static_cast<double>(9.0)));
        m_m = m_l;
        m_W_zero = W(Vector3D());
    }

public:
    /**
     * W(r,h) = (315/(64 pi h^9))(h^2-|r|^2)^3
     *        = (315/(64 pi h^9))(h^2-r*r)^3
     */
    static double W(const double r) {
        double res = 0.0;
        const double r2 = r * r;
        const double radius2 = m_radius * m_radius;
        if (r2 <= radius2) {
            res = pow(radius2 - r2, static_cast<double>(3.0)) * m_k;
        }
        return res;
    }

    static double W(const Vector3D &r) {
        double res = 0.0;
        const double r2 = r.lengthSquared();
        const double radius2 = m_radius * m_radius;
        if (r2 <= radius2) {
            res = pow(radius2 - r2, static_cast<double>(3.0)) * m_k;
        }
        return res;
    }

    /**
     * grad(W(r,h)) = r(-945/(32 pi h^9))(h^2-|r|^2)^2
     *              = r(-945/(32 pi h^9))(h^2-r*r)^2
     */
    static Vector3D gradW(const Vector3D &r) {
        Vector3D res;
        const double r2 = r.lengthSquared();
        const double radius2 = m_radius * m_radius;
        if (r2 <= radius2) {
            double tmp = radius2 - r2;
            res = m_l * tmp * tmp * r;
        } else
            res.setZero();

        return res;
    }

    /**
     * laplacian(W(r,h)) = (-945/(32 pi h^9))(h^2-|r|^2)(-7|r|^2+3h^2)
     *                   = (-945/(32 pi h^9))(h^2-r*r)(3 h^2-7 r*r)
     */
    static double laplacianW(const Vector3D &r) {
        double res;
        const double r2 = r.lengthSquared();
        const double radius2 = m_radius * m_radius;
        if (r2 <= radius2) {
            double tmp = radius2 - r2;
            double tmp2 = 3 * radius2 - 7 * r2;
            res = m_m * tmp * tmp2;
        } else
            res = 0.;

        return res;
    }

    static double W_zero() { return m_W_zero; }
};

// MARK: - SpikyKernel
/**
 * \brief Spiky kernel.
 */
class SpikyKernel {
protected:
    static double m_radius;
    static double m_k;
    static double m_l;
    static double m_W_zero;

public:
    static double getRadius() { return m_radius; }
    static void setRadius(double val) {
        m_radius = val;
        const double radius6 = pow(m_radius, static_cast<double>(6.0));
        const auto pi = static_cast<double>(M_PI);
        m_k = static_cast<double>(15.0) / (pi * radius6);
        m_l = -static_cast<double>(45.0) / (pi * radius6);
        m_W_zero = W(Vector3D());
    }

public:
    /**
     * W(r,h) = 15/(pi*h^6) * (h-r)^3
     */
    static double W(const double r) {
        double res = 0.0;
        const double r2 = r * r;
        const double radius2 = m_radius * m_radius;
        if (r2 <= radius2) {
            const double hr3 = pow(m_radius - r, static_cast<double>(3.0));
            res = m_k * hr3;
        }
        return res;
    }

    static double W(const Vector3D &r) {
        double res = 0.0;
        const double r2 = r.lengthSquared();
        const double radius2 = m_radius * m_radius;
        if (r2 <= radius2) {
            const double hr3 = pow(m_radius - sqrt(r2), static_cast<double>(3.0));
            res = m_k * hr3;
        }
        return res;
    }

    /**
     * grad(W(r,h)) = -r(45/(pi*h^6) * (h-r)^2)
     */
    static Vector3D gradW(const Vector3D &r) {
        Vector3D res;
        const double r2 = r.lengthSquared();
        const double radius2 = m_radius * m_radius;
        if (r2 <= radius2) {
            const double r_l = sqrt(r2);
            const double hr = m_radius - r_l;
            const double hr2 = hr * hr;
            res = m_l * hr2 * r * (static_cast<double>(1.0) / r_l);
        } else
            res.setZero();

        return res;
    }

    static double W_zero() { return m_W_zero; }
};

// MARK: - WendlandQuinticC2Kernel
/**
 * \brief quintic Wendland C2 kernel.
 */
class WendlandQuinticC2Kernel {
protected:
    static double m_radius;
    static double m_k;
    static double m_l;
    static double m_W_zero;

public:
    static double getRadius() { return m_radius; }
    static void setRadius(double val) {
        m_radius = val;
        const auto pi = static_cast<double>(M_PI);

        const double h3 = m_radius * m_radius * m_radius;
        m_k = static_cast<double>(21.0) / (static_cast<double>(2.0) * pi * h3);
        m_l = -static_cast<double>(210.0) / (pi * h3);
        m_W_zero = W(0.0);
    }

public:
    static double W(const double r) {
        double res = 0.0;
        const double q = r / m_radius;
        if (q <= 1.0)
            res = m_k * pow(static_cast<double>(1.0) - q, static_cast<double>(4.0)) *
                  (static_cast<double>(4.0) * q + static_cast<double>(1.0));
        return res;
    }

    static double W(const Vector3D &r) { return W(r.length()); }

    static Vector3D gradW(const Vector3D &r) {
        Vector3D res;
        const double rl = r.length();
        const double q = rl / m_radius;
        if (q <= 1.0) {
            const Vector3D gradq = r * (static_cast<double>(1.0) / (rl * m_radius));
            res = m_l * gradq * pow(static_cast<double>(1.0) - q, 3);
        } else
            res.setZero();

        return res;
    }

    static double W_zero() { return m_W_zero; }
};

// MARK: - CohesionKernel
/**
 * \brief Cohesion kernel used for the surface tension method of Akinci el al. [ATT13].
 *
 * References:
 * - [AAT13] Nadir Akinci, Gizem Akinci, and Matthias Teschner. Versatile surface tension and adhesion for sph fluids.
 * ACM Trans. Graph., 32(6):182:1-182:8, November 2013. URL: http://doi.acm.org/10.1145/2508363.2508395
 */
class CohesionKernel {
protected:
    static double m_radius;
    static double m_k;
    static double m_c;
    static double m_W_zero;

public:
    static double getRadius() { return m_radius; }
    static void setRadius(double val) {
        m_radius = val;
        const auto pi = static_cast<double>(M_PI);
        m_k = static_cast<double>(32.0) / (pi * pow(m_radius, static_cast<double>(9.0)));
        m_c = pow(m_radius, static_cast<double>(6.0)) / static_cast<double>(64.0);
        m_W_zero = W(Vector3D());
    }

public:
    /**
     * W(r,h) = (32/(pi h^9))(h-r)^3*r^3					if h/2 < r <= h
     *          (32/(pi h^9))(2*(h-r)^3*r^3 - h^6/64		if 0 < r <= h/2
     */
    static double W(const double r) {
        double res = 0.0;
        const double r2 = r * r;
        const double radius2 = m_radius * m_radius;
        if (r2 <= radius2) {
            const double r1 = sqrt(r2);
            const double r3 = r2 * r1;
            if (r1 > 0.5 * m_radius)
                res = m_k * pow(m_radius - r1, static_cast<double>(3.0)) * r3;
            else
                res = m_k * static_cast<double>(2.0) * pow(m_radius - r1, static_cast<double>(3.0)) * r3 - m_c;
        }
        return res;
    }

    static double W(const Vector3D &r) {
        double res = 0.0;
        const double r2 = r.lengthSquared();
        const double radius2 = m_radius * m_radius;
        if (r2 <= radius2) {
            const double r1 = sqrt(r2);
            const double r3 = r2 * r1;
            if (r1 > 0.5 * m_radius)
                res = m_k * pow(m_radius - r1, static_cast<double>(3.0)) * r3;
            else
                res = m_k * static_cast<double>(2.0) * pow(m_radius - r1, static_cast<double>(3.0)) * r3 - m_c;
        }
        return res;
    }

    static double W_zero() { return m_W_zero; }
};

// MARK: - AdhesionKernel
/**
 * \brief Adhesion kernel used for the surface tension method of Akinci el al. [ATT13].
 *
 * References:
 * - [AAT13] Nadir Akinci, Gizem Akinci, and Matthias Teschner. Versatile surface tension and adhesion for sph fluids.
 * ACM Trans. Graph., 32(6):182:1-182:8, November 2013. URL: http://doi.acm.org/10.1145/2508363.2508395
 */
class AdhesionKernel {
protected:
    static double m_radius;
    static double m_k;
    static double m_W_zero;

public:
    static double getRadius() { return m_radius; }
    static void setRadius(double val) {
        m_radius = val;
        m_k = static_cast<double>(0.007) / pow(m_radius, static_cast<double>(3.25));
        m_W_zero = W(Vector3D());
    }

public:
    /**
     * W(r,h) = (0.007/h^3.25)(-4r^2/h + 6r -2h)^0.25					if h/2 < r <= h
     */
    static double W(const double r) {
        double res = 0.0;
        const double r2 = r * r;
        const double radius2 = m_radius * m_radius;
        if (r2 <= radius2) {
            if (r > 0.5 * m_radius)
                res = m_k * pow(-static_cast<double>(4.0) * r2 / m_radius + static_cast<double>(6.0) * r -
                                        static_cast<double>(2.0) * m_radius,
                                static_cast<double>(0.25));
        }
        return res;
    }

    static double W(const Vector3D &r) {
        double res = 0.0;
        const double r2 = r.lengthSquared();
        const double radius2 = m_radius * m_radius;
        if (r2 <= radius2) {
            const double rl = sqrt(r2);
            if (rl > 0.5 * m_radius)
                res = m_k * pow(-static_cast<double>(4.0) * r2 / m_radius + static_cast<double>(6.0) * rl -
                                        static_cast<double>(2.0) * m_radius,
                                static_cast<double>(0.25));
        }
        return res;
    }

    static double W_zero() { return m_W_zero; }
};

// MARK: - CubicKernel2D
/**
 * \brief Cubic spline kernel (2D).
 */
class CubicKernel2D {
protected:
    static double m_radius;
    static double m_k;
    static double m_l;
    static double m_W_zero;

public:
    static double getRadius() { return m_radius; }
    static void setRadius(double val) {
        m_radius = val;
        const auto pi = static_cast<double>(M_PI);

        const double h2 = m_radius * m_radius;
        m_k = static_cast<double>(40.0) / (static_cast<double>(7.0) * (pi * h2));
        m_l = static_cast<double>(240.0) / (static_cast<double>(7.0) * (pi * h2));

        m_W_zero = W(Vector3D());
    }

public:
    static double W(const double r) {
        double res = 0.0;
        const double q = r / m_radius;
        if (q <= 1.0) {
            if (q <= 0.5) {
                const double q2 = q * q;
                const double q3 = q2 * q;
                res = m_k * (static_cast<double>(6.0) * q3 - static_cast<double>(6.0) * q2 + static_cast<double>(1.0));
            } else {
                res = m_k * (static_cast<double>(2.0) * pow(static_cast<double>(1.0) - q, static_cast<double>(3.0)));
            }
        }
        return res;
    }

    static double W(const Vector3D &r) { return W(r.length()); }

    static Vector3D gradW(const Vector3D &r) {
        Vector3D res;
        const double rl = r.length();
        const double q = rl / m_radius;
        if (q <= 1.0) {
            if (rl > 1.0e-5) {
                const Vector3D gradq = r * (static_cast<double>(1.0) / (rl * m_radius));
                if (q <= 0.5) {
                    res = m_l * q * (static_cast<double>(3.0) * q - static_cast<double>(2.0)) * gradq;
                } else {
                    const double factor = static_cast<double>(1.0) - q;
                    res = m_l * (-factor * factor) * gradq;
                }
            }
        } else
            res.setZero();

        return res;
    }

    static double W_zero() { return m_W_zero; }
};

// MARK: - WendlandQuinticC2Kernel2D
/**
 * \brief Wendland Quintic C2 spline kernel (2D).
 */
class WendlandQuinticC2Kernel2D {
protected:
    static double m_radius;
    static double m_k;
    static double m_l;
    static double m_W_zero;

public:
    static double getRadius() { return m_radius; }
    static void setRadius(double val) {
        m_radius = val;
        const auto pi = static_cast<double>(M_PI);

        const double h2 = m_radius * m_radius;
        m_k = static_cast<double>(7.0) / (pi * h2);
        m_l = -static_cast<double>(140.0) / (pi * h2);

        m_W_zero = W(Vector3D());
    }

public:
    static double W(const double r) {
        double res = 0.0;
        const double q = r / m_radius;
        if (q <= 1.0)
            res = m_k * pow(static_cast<double>(1.0) - q, static_cast<double>(4.0)) *
                  (static_cast<double>(4.0) * q + static_cast<double>(1.0));
        return res;
    }

    static double W(const Vector3D &r) { return W(r.length()); }

    static Vector3D gradW(const Vector3D &r) {
        Vector3D res;
        const double rl = r.length();
        const double q = rl / m_radius;
        if (q <= 1.0) {
            const Vector3D gradq = r * (static_cast<double>(1.0) / (rl * m_radius));
            res = m_l * q * pow(static_cast<double>(1.0) - q, static_cast<double>(3.0)) * gradq;
        } else
            res.setZero();

        return res;
    }

    static double W_zero() { return m_W_zero; }
};

// MARK: - PrecomputedKernel
/**
 * \brief Precomputed kernel which is based on a lookup table as described by Bender and Koschier [BK15,BK17].
 *
 * The lookup tables can be used in combination with any kernel.
 *
 * References:
 * - [BK15] Jan Bender and Dan Koschier. Divergence-free smoothed particle hydrodynamics. In ACM SIGGRAPH / Eurographics
 * Symposium on Computer Animation, SCA '15, 147-155. New York, NY, USA, 2015. ACM. URL:
 * http://doi.acm.org/10.1145/2786784.2786796
 * - [BK17] Jan Bender and Dan Koschier. Divergence-free SPH for incompressible and viscous fluids. IEEE Transactions on
 * Visualization and Computer Graphics, 23(3):1193-1206, 2017. URL: http://dx.doi.org/10.1109/TVCG.2016.2578335
 */
template <typename KernelType, unsigned int resolution = 10000u>
class PrecomputedKernel {
protected:
    static double m_W[resolution];
    static double m_gradW[resolution + 1];
    static double m_radius;
    static double m_radius2;
    static double m_invStepSize;
    static double m_W_zero;

public:
    static double getRadius() { return m_radius; }
    static void setRadius(double val) {
        m_radius = val;
        m_radius2 = m_radius * m_radius;
        KernelType::setRadius(val);
        const double stepSize = m_radius / (double)(resolution - 1);
        m_invStepSize = static_cast<double>(1.0) / stepSize;
        for (unsigned int i = 0; i < resolution; i++) {
            const double posX = stepSize * (double)i;  // Store kernel values in the middle of an interval
            m_W[i] = KernelType::W(posX);
            KernelType::setRadius(val);
            if (posX > 1.0e-9)
                m_gradW[i] = KernelType::gradW(Vector3D(posX, 0.0, 0.0))[0] / posX;
            else
                m_gradW[i] = 0.0;
        }
        m_gradW[resolution] = 0.0;
        m_W_zero = W(static_cast<double>(0));
    }

public:
    static double W(const Vector3D &r) {
        double res = 0.0;
        const double r2 = r.lengthSquared();
        if (r2 <= m_radius2) {
            const double rl = sqrt(r2);
            const unsigned int pos = std::min<unsigned int>((unsigned int)(rl * m_invStepSize), resolution - 2u);
            res = static_cast<double>(0.5) * (m_W[pos] + m_W[pos + 1]);
        }
        return res;
    }

    static double W(const double r) {
        double res = 0.0;
        if (r <= m_radius) {
            const unsigned int pos = std::min<unsigned int>((unsigned int)(r * m_invStepSize), resolution - 2u);
            res = static_cast<double>(0.5) * (m_W[pos] + m_W[pos + 1]);
        }
        return res;
    }

    static Vector3D gradW(const Vector3D &r) {
        Vector3D res;
        const double rl = r.length();
        if (rl <= m_radius) {
            // const double rl = sqrt(r2);
            const unsigned int pos =
                    std::min<unsigned int>(static_cast<unsigned int>(rl * m_invStepSize), resolution - 2u);
            res = static_cast<double>(0.5) * (m_gradW[pos] + m_gradW[pos + 1]) * r;
        } else
            res.setZero();

        return res;
    }

    static double W_zero() { return m_W_zero; }
};

template <typename KernelType, unsigned int resolution>
double PrecomputedKernel<KernelType, resolution>::m_radius;
template <typename KernelType, unsigned int resolution>
double PrecomputedKernel<KernelType, resolution>::m_radius2;
template <typename KernelType, unsigned int resolution>
double PrecomputedKernel<KernelType, resolution>::m_W[resolution];
template <typename KernelType, unsigned int resolution>
double PrecomputedKernel<KernelType, resolution>::m_gradW[resolution + 1];
template <typename KernelType, unsigned int resolution>
double PrecomputedKernel<KernelType, resolution>::m_invStepSize;
template <typename KernelType, unsigned int resolution>
double PrecomputedKernel<KernelType, resolution>::m_W_zero;

}  // namespace vox

#include "vox.force/sph_kernels3-inl.h"
