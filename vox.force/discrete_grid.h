//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <array>
#include <fstream>
#include <vector>

#include "vox.math/bounding_box3.h"
#include "vox.math/matrix.h"
#include "vox.math/size3.h"
#include "vox.math/vector3.h"

namespace vox::flex {

class DiscreteGrid {
public:
    using CoefficientVector = Matrix<double, 32, 1>;
    using ContinuousFunction = std::function<double(Point3D const &)>;
    using Predicate = std::function<bool(Point3D const &, double)>;
    using SamplePredicate = std::function<bool(Point3D const &)>;

    DiscreteGrid() = default;
    DiscreteGrid(BoundingBox3D const &domain, Size3 const &resolution)
        : m_domain(domain), m_resolution(resolution), m_n_fields(0u) {
        m_cell_size =
                Vector3D(domain.width() / resolution.x, domain.height() / resolution.y, domain.depth() / resolution.z);
        m_inv_cell_size = Vector3D(1 / m_cell_size.x, 1 / m_cell_size.y, 1 / m_cell_size.z);
        m_n_cells = resolution.x * resolution.y * resolution.z;
    }
    virtual ~DiscreteGrid() = default;

    virtual void save(std::string const &filename) const = 0;
    virtual void load(std::string const &filename) = 0;

    virtual unsigned int addFunction(ContinuousFunction const &func,
                                     bool verbose = false,
                                     SamplePredicate const &pred = nullptr) = 0;

    double interpolate(Point3D const &xi, Vector3D *gradient = nullptr) const { return interpolate(0u, xi, gradient); }

    virtual double interpolate(unsigned int field_id, Point3D const &xi, Vector3D *gradient = nullptr) const = 0;

    /**
     * @brief Determines the shape functions for the discretization with ID
     * field_id at point xi.
     *
     * @param field_id Discretization ID
     * @param x Location where the shape functions should be determined
     * @param cell cell of x
     * @param c0 vector required for the interpolation
     * @param N	shape functions for the cell of x
     * @param dN (Optional) derivatives of the shape functions, required to
     * compute the gradient
     * @return Success of the function.
     */
    virtual bool determineShapeFunctions(unsigned int field_id,
                                         Point3D const &x,
                                         std::array<unsigned int, 32> &cell,
                                         Vector3D &c0,
                                         Matrix<double, 32, 1> &N,
                                         Matrix<double, 32, 3> *dN = nullptr) const = 0;

    /**
     * @brief Evaluates the given discretization with ID field_id at point xi.
     *
     * @param field_id Discretization ID
     * @param xi Location where the discrete function is evaluated
     * @param cell cell of xi
     * @param c0 vector required for the interpolation
     * @param N	shape functions for the cell of xi
     * @param gradient (Optional) if a pointer to a vector is passed the gradient
     * of the discrete function will be evaluated
     * @param dN (Optional) derivatives of the shape functions, required to
     * compute the gradient
     * @return double Results of the evaluation of the discrete function at point
     * xi
     */
    virtual double interpolate(unsigned int field_id,
                               Point3D const &xi,
                               const std::array<unsigned int, 32> &cell,
                               const Vector3D &c0,
                               const Matrix<double, 32, 1> &N,
                               Vector3D *gradient = nullptr,
                               Matrix<double, 32, 3> *dN = nullptr) const = 0;

    virtual void reduceField(unsigned int field_id, Predicate pred) {}

    [[nodiscard]] Size3 singleToMultiIndex(unsigned int i) const;
    [[nodiscard]] size_t multiToSingleIndex(Size3 const &ijk) const;

    [[nodiscard]] BoundingBox3D subdomain(Size3 const &ijk) const;
    [[nodiscard]] BoundingBox3D subdomain(unsigned int l) const;

    [[nodiscard]] BoundingBox3D const &domain() const { return m_domain; }
    [[nodiscard]] Size3 const &resolution() const { return m_resolution; };
    [[nodiscard]] Vector3D const &cellSize() const { return m_cell_size; }
    [[nodiscard]] Vector3D const &invCellSize() const { return m_inv_cell_size; }

protected:
    BoundingBox3D m_domain;
    Size3 m_resolution;
    Vector3D m_cell_size;
    Vector3D m_inv_cell_size;
    std::size_t m_n_cells{};
    std::size_t m_n_fields{};
};
}  // namespace vox::flex
