//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/sph_kernels3.h"

namespace vox {
double CubicKernel::m_radius;
double CubicKernel::m_k;
double CubicKernel::m_l;
double CubicKernel::m_W_zero;

double Poly6Kernel::m_radius;
double Poly6Kernel::m_k;
double Poly6Kernel::m_l;
double Poly6Kernel::m_m;
double Poly6Kernel::m_W_zero;

double SpikyKernel::m_radius;
double SpikyKernel::m_k;
double SpikyKernel::m_l;
double SpikyKernel::m_W_zero;

double WendlandQuinticC2Kernel::m_radius;
double WendlandQuinticC2Kernel::m_k;
double WendlandQuinticC2Kernel::m_l;
double WendlandQuinticC2Kernel::m_W_zero;

double CohesionKernel::m_radius;
double CohesionKernel::m_k;
double CohesionKernel::m_c;
double CohesionKernel::m_W_zero;

double AdhesionKernel::m_radius;
double AdhesionKernel::m_k;
double AdhesionKernel::m_W_zero;

double CubicKernel2D::m_radius;
double CubicKernel2D::m_k;
double CubicKernel2D::m_l;
double CubicKernel2D::m_W_zero;

double WendlandQuinticC2Kernel2D::m_radius;
double WendlandQuinticC2Kernel2D::m_k;
double WendlandQuinticC2Kernel2D::m_l;
double WendlandQuinticC2Kernel2D::m_W_zero;

}  // namespace vox
