//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <vector>

#include "vox.geometry/array2.h"
#include "vox.geometry/array_samplers2.h"
#include "vox.geometry/vector_grid2.h"

namespace vox {

//! \brief Abstract base class for 2-D collocated vector grid structure.
class CollocatedVectorGrid2 : public VectorGrid2 {
public:
    //! Constructs an empty grid.
    CollocatedVectorGrid2();

    //! Default destructor.
    ~CollocatedVectorGrid2() override;

    //! Returns the actual data point size.
    [[nodiscard]] virtual Size2 dataSize() const = 0;

    //!
    //! \brief Returns data position for the grid point at (0, 0).
    //!
    //! Note that this is different from origin() since origin() returns
    //! the lower corner point of the bounding box.
    //!
    [[nodiscard]] virtual Point2D dataOrigin() const = 0;

    //! Returns the grid data at given data point.
    const Vector2D &operator()(size_t i, size_t j) const;

    //! Returns the grid data at given data point.
    Vector2D &operator()(size_t i, size_t j);

    //! Returns divergence at data point location.
    [[nodiscard]] double divergenceAtDataPoint(size_t i, size_t j) const;

    //! Returns curl at data point location.
    [[nodiscard]] double curlAtDataPoint(size_t i, size_t j) const;

    //! Returns the read-write data array accessor.
    VectorDataAccessor dataAccessor();

    //! Returns the read-only data array accessor.
    [[nodiscard]] ConstVectorDataAccessor constDataAccessor() const;

    //! Returns the function that maps data point to its position.
    [[nodiscard]] DataPositionFunc dataPosition() const;

    //!
    //! \brief Invokes the given function \p func for each data point.
    //!
    //! This function invokes the given function object \p func for each data
    //! point in serial manner. The input parameters are i and j indices of a
    //! data point. The order of execution is i-first, j-last.
    //!
    void forEachDataPointIndex(const std::function<void(size_t, size_t)> &func) const;

    //!
    //! \brief Invokes the given function \p func for each data point
    //! parallelly.
    //!
    //! This function invokes the given function object \p func for each data
    //! point in parallel manner. The input parameters are i and j indices of a
    //! data point. The order of execution can be arbitrary since it's
    //! multi-threaded.
    //!
    void parallelForEachDataPointIndex(const std::function<void(size_t, size_t)> &func) const;

    // VectorField2 implementations

    //! Returns sampled value at given position \p x.
    [[nodiscard]] Vector2D sample(const Point2D &x) const override;

    //! Returns divergence at given position \p x.
    [[nodiscard]] double divergence(const Point2D &x) const override;

    //! Returns curl at given position \p x.
    [[nodiscard]] double curl(const Point2D &x) const override;

    //!
    //! \brief Returns the sampler function.
    //!
    //! This function returns the data sampler function object. The sampling
    //! function is linear.
    //!
    [[nodiscard]] std::function<Vector2D(const Point2D &)> sampler() const override;

protected:
    //! Swaps the data storage and predefined samplers with given grid.
    void swapCollocatedVectorGrid(CollocatedVectorGrid2 *other);

    //! Sets the data storage and predefined samplers with given grid.
    void setCollocatedVectorGrid(const CollocatedVectorGrid2 &other);

    //! Fetches the data into a continuous linear array.
    void getData(std::vector<double> *data) const override;

    //! Sets the data from a continuous linear array.
    void setData(const std::vector<double> &data) override;

private:
    Array2<Vector2D> _data;
    LinearArraySampler2<Vector2D, double> _linearSampler;
    std::function<Vector2D(const Point2D &)> _sampler;

    void onResize(const Size2 &resolution,
                  const Vector2D &gridSpacing,
                  const Point2D &origin,
                  const Vector2D &initialValue) final;

    void resetSampler();
};

//! Shared pointer for the CollocatedVectorGrid2 type.
typedef std::shared_ptr<CollocatedVectorGrid2> CollocatedVectorGrid2Ptr;

}  // namespace vox