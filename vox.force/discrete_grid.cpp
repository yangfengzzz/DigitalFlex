//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.force/discrete_grid.h"

namespace vox::flex {

Size3 DiscreteGrid::singleToMultiIndex(unsigned int l) const {
    auto n01 = m_resolution.x * m_resolution.y;
    auto k = l / n01;
    auto temp = l % n01;
    auto j = temp / m_resolution.x;
    auto i = temp % m_resolution.x;
    return {i, j, k};
}

size_t DiscreteGrid::multiToSingleIndex(Size3 const &ijk) const {
    return m_resolution.y * m_resolution.x * ijk.z + m_resolution.x * ijk.y + ijk.x;
}

BoundingBox3D DiscreteGrid::subdomain(Size3 const &ijk) const {
    auto origin = m_domain.lower_corner + Vector3D(ijk.x * m_cell_size.x, ijk.y * m_cell_size.y, ijk.z * m_cell_size.z);
    return {origin, origin + m_cell_size};
}

BoundingBox3D DiscreteGrid::subdomain(unsigned int l) const { return subdomain(singleToMultiIndex(l)); }

}  // namespace vox::flex
