// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.geometry/vertex_centered_vector_grid2.h"

#include <utility>

#include "vox.geometry/array_samplers2.h"
#include "vox.geometry/parallel.h"
#include "vox.geometry/private_helpers.h"

using namespace vox;

VertexCenteredVectorGrid2::VertexCenteredVectorGrid2() = default;

VertexCenteredVectorGrid2::VertexCenteredVectorGrid2(size_t resolutionX,
                                                     size_t resolutionY,
                                                     double gridSpacingX,
                                                     double gridSpacingY,
                                                     double originX,
                                                     double originY,
                                                     double initialValueU,
                                                     double initialValueV) {
    resize(resolutionX, resolutionY, gridSpacingX, gridSpacingY, originX, originY, initialValueU, initialValueV);
}

VertexCenteredVectorGrid2::VertexCenteredVectorGrid2(const Size2 &resolution,
                                                     const Vector2D &gridSpacing,
                                                     const Point2D &origin,
                                                     const Vector2D &initialValue) {
    resize(resolution, gridSpacing, origin, initialValue);
}

VertexCenteredVectorGrid2::VertexCenteredVectorGrid2(const VertexCenteredVectorGrid2 &other) { set(other); }

Size2 VertexCenteredVectorGrid2::dataSize() const {
    if (resolution() != Size2(0, 0)) {
        return resolution() + Size2(1, 1);
    } else {
        return {0, 0};
    }
}

Point2D VertexCenteredVectorGrid2::dataOrigin() const { return origin(); }

void VertexCenteredVectorGrid2::swap(Grid2 *other) {
    auto *sameType = dynamic_cast<VertexCenteredVectorGrid2 *>(other);
    if (sameType != nullptr) {
        swapCollocatedVectorGrid(sameType);
    }
}

void VertexCenteredVectorGrid2::set(const VertexCenteredVectorGrid2 &other) { setCollocatedVectorGrid(other); }

VertexCenteredVectorGrid2 &VertexCenteredVectorGrid2::operator=(const VertexCenteredVectorGrid2 &other) {
    set(other);
    return *this;
}

void VertexCenteredVectorGrid2::fill(const Vector2D &value, ExecutionPolicy policy) {
    Size2 size = dataSize();
    auto acc = dataAccessor();
    parallelFor(
            kZeroSize, size.x, kZeroSize, size.y, [value, &acc](size_t i, size_t j) { acc(i, j) = value; }, policy);
}

void VertexCenteredVectorGrid2::fill(const std::function<Vector2D(const Point2D &)> &func, ExecutionPolicy policy) {
    Size2 size = dataSize();
    auto acc = dataAccessor();
    DataPositionFunc pos = dataPosition();
    parallelFor(
            kZeroSize, size.x, kZeroSize, size.y,
            [&func, &acc, &pos](size_t i, size_t j) { acc(i, j) = func(pos(i, j)); }, policy);
}

std::shared_ptr<VectorGrid2> VertexCenteredVectorGrid2::clone() const {
        return CLONE_W_CUSTOM_DELETER(VertexCenteredVectorGrid2)}

VertexCenteredVectorGrid2::Builder VertexCenteredVectorGrid2::builder() {
    return {};
}

VertexCenteredVectorGrid2::Builder &VertexCenteredVectorGrid2::Builder::withResolution(const Size2 &resolution) {
    _resolution = resolution;
    return *this;
}

VertexCenteredVectorGrid2::Builder &VertexCenteredVectorGrid2::Builder::withResolution(size_t resolutionX,
                                                                                       size_t resolutionY) {
    _resolution.x = resolutionX;
    _resolution.y = resolutionY;
    return *this;
}

VertexCenteredVectorGrid2::Builder &VertexCenteredVectorGrid2::Builder::withGridSpacing(const Vector2D &gridSpacing) {
    _gridSpacing = gridSpacing;
    return *this;
}

VertexCenteredVectorGrid2::Builder &VertexCenteredVectorGrid2::Builder::withGridSpacing(double gridSpacingX,
                                                                                        double gridSpacingY) {
    _gridSpacing.x = gridSpacingX;
    _gridSpacing.y = gridSpacingY;
    return *this;
}

VertexCenteredVectorGrid2::Builder &VertexCenteredVectorGrid2::Builder::withOrigin(const Point2D &gridOrigin) {
    _gridOrigin = gridOrigin;
    return *this;
}

VertexCenteredVectorGrid2::Builder &VertexCenteredVectorGrid2::Builder::withOrigin(double gridOriginX,
                                                                                   double gridOriginY) {
    _gridOrigin.x = gridOriginX;
    _gridOrigin.y = gridOriginY;
    return *this;
}

VertexCenteredVectorGrid2::Builder &VertexCenteredVectorGrid2::Builder::withInitialValue(const Vector2D &initialVal) {
    _initialVal = initialVal;
    return *this;
}

VertexCenteredVectorGrid2::Builder &VertexCenteredVectorGrid2::Builder::withInitialValue(double initialValX,
                                                                                         double initialValY) {
    _initialVal.x = initialValX;
    _initialVal.y = initialValY;
    return *this;
}

VertexCenteredVectorGrid2 VertexCenteredVectorGrid2::Builder::build() const {
    return VertexCenteredVectorGrid2(_resolution, _gridSpacing, _gridOrigin, _initialVal);
}

VertexCenteredVectorGrid2Ptr VertexCenteredVectorGrid2::Builder::makeShared() const {
    return {new VertexCenteredVectorGrid2(_resolution, _gridSpacing, _gridOrigin, _initialVal),
            [](VertexCenteredVectorGrid2 *obj) { delete obj; }};
}

VectorGrid2Ptr VertexCenteredVectorGrid2::Builder::build(const Size2 &resolution,
                                                         const Vector2D &gridSpacing,
                                                         const Point2D &gridOrigin,
                                                         const Vector2D &initialVal) const {
    return std::shared_ptr<VertexCenteredVectorGrid2>(
            new VertexCenteredVectorGrid2(resolution, gridSpacing, gridOrigin, initialVal),
            [](VertexCenteredVectorGrid2 *obj) { delete obj; });
}
