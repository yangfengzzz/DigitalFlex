// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "manual_tests.h"
#include "vox.force/sph_system_data3.h"
#include "vox.geometry/array2.h"
#include "vox.geometry/bcc_lattice_point_generator.h"
#include "vox.geometry/point_hash_grid_searcher2.h"
#include "vox.geometry/point_hash_grid_searcher3.h"
#include "vox.geometry/point_parallel_hash_grid_searcher2.h"
#include "vox.geometry/point_parallel_hash_grid_searcher3.h"
#include "vox.geometry/triangle_point_generator.h"
#include "vox.math/bounding_box2.h"
#include "vox.math/bounding_box3.h"

using namespace vox;

JET_TESTS(PointHashGridSearcher2);

JET_BEGIN_TEST_F(PointHashGridSearcher2, Build) {
    Array1<Point2D> points;
    TrianglePointGenerator pointsGenerator;
    BoundingBox2D bbox(Point2D(0, 0), Point2D(1, 1));
    double spacing = 0.1;

    pointsGenerator.generate(bbox, spacing, &points);

    PointHashGridSearcher2 pointSearcher(4, 4, 0.18);
    pointSearcher.build(ArrayAccessor1<Point2D>(points.size(), points.data()));

    Array2<double> grid(4, 4, 0.0);

    for (size_t j = 0; j < grid.size().y; ++j) {
        for (size_t i = 0; i < grid.size().x; ++i) {
            size_t key =
                    pointSearcher.getHashKeyFromBucketIndex(Point2I(static_cast<ssize_t>(i), static_cast<ssize_t>(j)));
            size_t value = pointSearcher.buckets()[key].size();
            grid(i, j) += static_cast<double>(value);
        }
    }

    saveData(grid.constAccessor(), "data_#grid2.npy");
}
JET_END_TEST_F

JET_TESTS(PointHashGridSearcher3);

JET_BEGIN_TEST_F(PointHashGridSearcher3, Build) {
    Array1<Point3D> points;
    BccLatticePointGenerator pointsGenerator;
    BoundingBox3D bbox(Point3D(0, 0, 0), Point3D(1, 1, 1));
    double spacing = 0.1;

    pointsGenerator.generate(bbox, spacing, &points);

    PointHashGridSearcher3 pointSearcher(4, 4, 4, 0.18);
    pointSearcher.build(ArrayAccessor1<Point3D>(points.size(), points.data()));

    Array2<double> grid(4, 4, 0.0);

    for (size_t j = 0; j < grid.size().y; ++j) {
        for (size_t i = 0; i < grid.size().x; ++i) {
            size_t key = pointSearcher.getHashKeyFromBucketIndex(
                    Point3I(static_cast<ssize_t>(i), static_cast<ssize_t>(j), 0));
            size_t value = pointSearcher.buckets()[key].size();
            grid(i, j) += static_cast<double>(value);
        }
    }

    saveData(grid.constAccessor(), "data_#grid2.npy");
}
JET_END_TEST_F

JET_TESTS(PointParallelHashGridSearcher2);

JET_BEGIN_TEST_F(PointParallelHashGridSearcher2, Build) {
    Array1<Point2D> points;
    TrianglePointGenerator pointsGenerator;
    BoundingBox2D bbox(Point2D(0, 0), Point2D(1, 1));
    double spacing = 0.1;

    pointsGenerator.generate(bbox, spacing, &points);

    PointParallelHashGridSearcher2 pointSearcher(4, 4, 0.18);
    pointSearcher.build(ArrayAccessor1<Point2D>(points.size(), points.data()));

    Array2<double> grid(4, 4, 0.0);

    for (size_t j = 0; j < grid.size().y; ++j) {
        for (size_t i = 0; i < grid.size().x; ++i) {
            size_t key =
                    pointSearcher.getHashKeyFromBucketIndex(Point2I(static_cast<ssize_t>(i), static_cast<ssize_t>(j)));
            size_t start = pointSearcher.startIndexTable()[key];
            size_t end = pointSearcher.endIndexTable()[key];
            size_t value = end - start;
            grid(i, j) += static_cast<double>(value);
        }
    }

    saveData(grid.constAccessor(), "data_#grid2.npy");
}
JET_END_TEST_F

JET_TESTS(PointParallelHashGridSearcher3);

JET_BEGIN_TEST_F(PointParallelHashGridSearcher3, Build) {
    Array1<Point3D> points;
    BccLatticePointGenerator pointsGenerator;
    BoundingBox3D bbox(Point3D(0, 0, 0), Point3D(1, 1, 1));
    double spacing = 0.1;

    pointsGenerator.generate(bbox, spacing, &points);

    PointParallelHashGridSearcher3 pointSearcher(4, 4, 4, 0.18);
    pointSearcher.build(ArrayAccessor1<Point3D>(points.size(), points.data()));

    Array2<double> grid(4, 4, 0.0);

    for (size_t j = 0; j < grid.size().y; ++j) {
        for (size_t i = 0; i < grid.size().x; ++i) {
            size_t key = pointSearcher.getHashKeyFromBucketIndex(
                    Point3I(static_cast<ssize_t>(i), static_cast<ssize_t>(j), 0));
            size_t start = pointSearcher.startIndexTable()[key];
            size_t end = pointSearcher.endIndexTable()[key];
            size_t value = end - start;
            grid(i, j) += static_cast<double>(value);
        }
    }

    saveData(grid.constAccessor(), "data_#grid2.npy");
}
JET_END_TEST_F
