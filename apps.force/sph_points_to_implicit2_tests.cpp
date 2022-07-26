// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include <random>

#include "manual_tests.h"
#include "vox.force/sph_points_to_implicit2.h"
#include "vox.geometry/cell_centered_scalar_grid2.h"

using namespace vox;

JET_TESTS(SphPointsToImplicit2);

JET_BEGIN_TEST_F(SphPointsToImplicit2, ConvertTwo) {
    Array1<Point2D> points;

    std::mt19937 rng{0};
    std::uniform_real_distribution<> dist(0.2, 0.8);
    for (size_t i = 0; i < 2; ++i) {
        points.append({dist(rng), dist(rng)});
    }

    CellCenteredScalarGrid2 grid(512, 512, 1.0 / 512, 1.0 / 512);

    SphPointsToImplicit2 converter(0.1);
    converter.convert(points.constAccessor(), &grid);

    saveData(grid.constDataAccessor(), "data_#grid2.npy");
    saveData(grid.constDataAccessor(), "data_#grid2,iso.npy");
}
JET_END_TEST_F

JET_BEGIN_TEST_F(SphPointsToImplicit2, ConvertMany) {
    Array1<Point2D> points;

    std::mt19937 rng{0};
    std::uniform_real_distribution<> dist(0.2, 0.8);
    for (size_t i = 0; i < 200; ++i) {
        points.append({dist(rng), dist(rng)});
    }

    CellCenteredScalarGrid2 grid(512, 512, 1.0 / 512, 1.0 / 512);

    SphPointsToImplicit2 converter(0.1);
    converter.convert(points.constAccessor(), &grid);

    saveData(grid.constDataAccessor(), "data_#grid2.npy");
    saveData(grid.constDataAccessor(), "data_#grid2,iso.npy");
}
JET_END_TEST_F
