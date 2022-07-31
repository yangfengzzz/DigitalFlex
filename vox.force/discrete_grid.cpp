//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/discrete_grid.h"

namespace vox::flex {

DiscreteGrid::MultiIndex DiscreteGrid::singleToMultiIndex(unsigned int l) const {
    auto n01 = m_resolution[0] * m_resolution[1];
    auto k = l / n01;
    auto temp = l % n01;
    auto j = temp / m_resolution[0];
    auto i = temp % m_resolution[0];
    return {{i, j, k}};
}

size_t DiscreteGrid::multiToSingleIndex(MultiIndex const &ijk) const {
    return m_resolution[1] * m_resolution[0] * ijk[2] + m_resolution[0] * ijk[1] + ijk[0];
}

BoundingBox3D DiscreteGrid::subdomain(MultiIndex const &ijk) const {
    auto origin =
            m_domain.lower_corner + Vector3D(ijk[0] * m_cell_size.x, ijk[1] * m_cell_size.y, ijk[2] * m_cell_size.z);
    return {origin, origin + m_cell_size};
}

BoundingBox3D DiscreteGrid::subdomain(unsigned int l) const { return subdomain(singleToMultiIndex(l)); }

}  // namespace vox::flex
