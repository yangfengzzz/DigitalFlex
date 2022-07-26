// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include <random>

#include "manual_tests.h"
#include "vox.base/parallel.h"
#include "vox.force/sph_system_data2.h"
#include "vox.force/sph_system_data3.h"
#include "vox.geometry/bcc_lattice_point_generator.h"
#include "vox.geometry/cell_centered_scalar_grid2.h"
#include "vox.geometry/triangle_point_generator.h"
#include "vox.math/bounding_box3.h"

using namespace vox;

JET_TESTS(SphSystemData2);

JET_BEGIN_TEST_F(SphSystemData2, Interpolate) {
    Array1<Point2D> points;
    TrianglePointGenerator pointsGenerator;
    BoundingBox2D bbox(Point2D(0, 0), Point2D(1, 1));
    double spacing = 0.1;

    pointsGenerator.generate(bbox, spacing, &points);

    SphSystemData2 sphSystem;
    sphSystem.addParticles(ConstArrayAccessor1<Point2D>(points.size(), points.data()));
    sphSystem.setTargetSpacing(spacing);
    sphSystem.buildNeighborSearcher();
    sphSystem.buildNeighborLists();
    sphSystem.updateDensities();

    Array1<double> data(points.size(), 1.0);

    CellCenteredScalarGrid2 grid(512, 512, 1.0 / 512, 1.0 / 512);

    auto gridPos = grid.dataPosition();
    parallelFor(kZeroSize, grid.dataSize().x, kZeroSize, grid.dataSize().y, [&](size_t i, size_t j) {
        Point2D xy = gridPos(i, j);
        Point2D p(xy.x, xy.y);
        grid(i, j) = sphSystem.interpolate(p, data.constAccessor());
    });

    saveData(grid.constDataAccessor(), "data_#grid2.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(SphSystemData2, Gradient) {
    Array1<Point2D> points;
    TrianglePointGenerator pointsGenerator;
    BoundingBox2D bbox(Point2D(0, 0), Point2D(1, 1));
    double spacing = 0.1;

    pointsGenerator.generate(bbox, spacing, &points);

    SphSystemData2 sphSystem;
    sphSystem.addParticles(ConstArrayAccessor1<Point2D>(points.size(), points.data()));
    sphSystem.setTargetSpacing(spacing);
    sphSystem.buildNeighborSearcher();
    sphSystem.buildNeighborLists();
    sphSystem.updateDensities();

    Array1<double> data(points.size()), gradX(points.size()), gradY(points.size());
    std::mt19937 rng(0);
    std::uniform_real_distribution<> d(0.0, 1.0);

    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = d(rng);
    }

    for (size_t i = 0; i < data.size(); ++i) {
        Vector2D g = sphSystem.gradientAt(i, data.constAccessor());
        gradX[i] = g.x;
        gradY[i] = g.y;
    }

    CellCenteredScalarGrid2 grid(64, 64, 1.0 / 64, 1.0 / 64);
    CellCenteredScalarGrid2 grid2(64, 64, 1.0 / 64, 1.0 / 64);

    auto gridPos = grid.dataPosition();
    parallelFor(kZeroSize, grid.dataSize().x, kZeroSize, grid.dataSize().y, [&](size_t i, size_t j) {
        Point2D xy = gridPos(i, j);
        Point2D p(xy.x, xy.y);
        grid(i, j) = sphSystem.interpolate(p, data.constAccessor());
    });

    saveData(grid.constDataAccessor(), "data_#grid2.npy");

    parallelFor(kZeroSize, grid.dataSize().x, kZeroSize, grid.dataSize().y, [&](size_t i, size_t j) {
        Point2D xy = gridPos(i, j);
        Point2D p(xy.x, xy.y);
        grid(i, j) = sphSystem.interpolate(p, gradX.constAccessor());
        grid2(i, j) = sphSystem.interpolate(p, gradY.constAccessor());
    });

    saveData(grid.constDataAccessor(), "gradient_#grid2,x.npy");
    saveData(grid2.constDataAccessor(), "gradient_#grid2,y.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(SphSystemData2, Laplacian) {
    Array1<Point2D> points;
    TrianglePointGenerator pointsGenerator;
    BoundingBox2D bbox(Point2D(0, 0), Point2D(1, 1));
    double spacing = 0.1;

    pointsGenerator.generate(bbox, spacing, &points);

    SphSystemData2 sphSystem;
    sphSystem.addParticles(ConstArrayAccessor1<Point2D>(points.size(), points.data()));
    sphSystem.setTargetSpacing(spacing);
    sphSystem.buildNeighborSearcher();
    sphSystem.buildNeighborLists();
    sphSystem.updateDensities();

    Array1<double> data(points.size()), laplacian(points.size());
    std::mt19937 rng(0);
    std::uniform_real_distribution<> d(0.0, 1.0);

    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = d(rng);
    }

    for (size_t i = 0; i < data.size(); ++i) {
        laplacian[i] = sphSystem.laplacianAt(i, data.constAccessor());
    }

    CellCenteredScalarGrid2 grid(512, 512, 1.0 / 512, 1.0 / 512);

    auto gridPos = grid.dataPosition();
    parallelFor(kZeroSize, grid.dataSize().x, kZeroSize, grid.dataSize().y, [&](size_t i, size_t j) {
        Point2D xy = gridPos(i, j);
        Point2D p(xy.x, xy.y);
        grid(i, j) = sphSystem.interpolate(p, data.constAccessor());
    });

    saveData(grid.constDataAccessor(), "data_#grid2.npy");

    parallelFor(kZeroSize, grid.dataSize().x, kZeroSize, grid.dataSize().y, [&](size_t i, size_t j) {
        Point2D xy = gridPos(i, j);
        Point2D p(xy.x, xy.y);
        grid(i, j) = sphSystem.interpolate(p, laplacian.constAccessor());
    });

    saveData(grid.constDataAccessor(), "laplacian_#grid2.npy");
}
JET_END_TEST_F

JET_TESTS(SphSystemData3);

JET_BEGIN_TEST_F(SphSystemData3, Interpolate) {
    Array1<Point3D> points;
    BccLatticePointGenerator pointsGenerator;
    BoundingBox3D bbox(Point3D(0, 0, 0), Point3D(1, 1, 1));
    double spacing = 0.1;

    pointsGenerator.generate(bbox, spacing, &points);

    SphSystemData3 sphSystem;
    sphSystem.addParticles(ConstArrayAccessor1<Point3D>(points.size(), points.data()));
    sphSystem.setTargetSpacing(spacing);
    sphSystem.buildNeighborSearcher();
    sphSystem.buildNeighborLists();
    sphSystem.updateDensities();

    Array1<double> data(points.size(), 1.0);

    CellCenteredScalarGrid2 grid(512, 512, 1.0 / 512, 1.0 / 512);

    auto gridPos = grid.dataPosition();
    parallelFor(kZeroSize, grid.dataSize().x, kZeroSize, grid.dataSize().y, [&](size_t i, size_t j) {
        Point2D xy = gridPos(i, j);
        Point3D p(xy.x, xy.y, 0.5);
        grid(i, j) = sphSystem.interpolate(p, data.constAccessor());
    });

    saveData(grid.constDataAccessor(), "data_#grid2.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(SphSystemData3, Gradient) {
    Array1<Point3D> points;
    BccLatticePointGenerator pointsGenerator;
    BoundingBox3D bbox(Point3D(0, 0, 0), Point3D(1, 1, 1));
    double spacing = 0.1;

    pointsGenerator.generate(bbox, spacing, &points);

    SphSystemData3 sphSystem;
    sphSystem.addParticles(ConstArrayAccessor1<Point3D>(points.size(), points.data()));
    sphSystem.setTargetSpacing(spacing);
    sphSystem.buildNeighborSearcher();
    sphSystem.buildNeighborLists();
    sphSystem.updateDensities();

    Array1<double> data(points.size());
    Array1<double> gradX(points.size()), gradY(points.size());
    std::mt19937 rng(0);
    std::uniform_real_distribution<> d(0.0, 1.0);

    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = d(rng);
    }

    for (size_t i = 0; i < data.size(); ++i) {
        Vector3D g = sphSystem.gradientAt(i, data.constAccessor());
        gradX[i] = g.x;
        gradY[i] = g.y;
    }

    CellCenteredScalarGrid2 grid(64, 64, 1.0 / 64, 1.0 / 64);
    CellCenteredScalarGrid2 grid2(64, 64, 1.0 / 64, 1.0 / 64);

    auto gridPos = grid.dataPosition();
    parallelFor(kZeroSize, grid.dataSize().x, kZeroSize, grid.dataSize().y, [&](size_t i, size_t j) {
        Point2D xy = gridPos(i, j);
        Point3D p(xy.x, xy.y, 0.5);
        grid(i, j) = sphSystem.interpolate(p, data.constAccessor());
    });

    saveData(grid.constDataAccessor(), "data_#grid2.npy");

    parallelFor(kZeroSize, grid.dataSize().x, kZeroSize, grid.dataSize().y, [&](size_t i, size_t j) {
        Point2D xy = gridPos(i, j);
        Point3D p(xy.x, xy.y, 0.5);
        grid(i, j) = sphSystem.interpolate(p, gradX.constAccessor());
        grid2(i, j) = sphSystem.interpolate(p, gradY.constAccessor());
    });

    saveData(grid.constDataAccessor(), "gradient_#grid2,x.npy");
    saveData(grid2.constDataAccessor(), "gradient_#grid2,y.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(SphSystemData3, Laplacian) {
    Array1<Point3D> points;
    BccLatticePointGenerator pointsGenerator;
    BoundingBox3D bbox(Point3D(0, 0, 0), Point3D(1, 1, 1));
    double spacing = 0.1;

    pointsGenerator.generate(bbox, spacing, &points);

    SphSystemData3 sphSystem;
    sphSystem.addParticles(ConstArrayAccessor1<Point3D>(points.size(), points.data()));
    sphSystem.setTargetSpacing(spacing);
    sphSystem.buildNeighborSearcher();
    sphSystem.buildNeighborLists();
    sphSystem.updateDensities();

    Array1<double> data(points.size()), laplacian(points.size());
    std::mt19937 rng(0);
    std::uniform_real_distribution<> d(0.0, 1.0);

    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = d(rng);
    }

    for (size_t i = 0; i < data.size(); ++i) {
        laplacian[i] = sphSystem.laplacianAt(i, data.constAccessor());
    }

    CellCenteredScalarGrid2 grid(512, 512, 1.0 / 512, 1.0 / 512);

    auto gridPos = grid.dataPosition();
    parallelFor(kZeroSize, grid.dataSize().x, kZeroSize, grid.dataSize().y, [&](size_t i, size_t j) {
        Point2D xy = gridPos(i, j);
        Point3D p(xy.x, xy.y, 0.5);
        grid(i, j) = sphSystem.interpolate(p, data.constAccessor());
    });

    saveData(grid.constDataAccessor(), "data_#grid2.npy");

    parallelFor(kZeroSize, grid.dataSize().x, kZeroSize, grid.dataSize().y, [&](size_t i, size_t j) {
        Point2D xy = gridPos(i, j);
        Point3D p(xy.x, xy.y, 0.5);
        grid(i, j) = sphSystem.interpolate(p, laplacian.constAccessor());
    });

    saveData(grid.constDataAccessor(), "laplacian_#grid2.npy");
}
JET_END_TEST_F
