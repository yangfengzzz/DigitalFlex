// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include <gtest/gtest.h>

#include "vox.force/pci_sph_solver2.h"

using namespace vox;

TEST(PciSphSolver2, UpdateEmpty) {
    // Empty solver test
    PciSphSolver2 solver;
    Frame frame(0, 0.01);
    solver.update(frame++);
    solver.update(frame);
}

TEST(PciSphSolver2, Parameters) {
    PciSphSolver2 solver;

    solver.setMaxDensityErrorRatio(5.0);
    EXPECT_DOUBLE_EQ(5.0, solver.maxDensityErrorRatio());

    solver.setMaxDensityErrorRatio(-1.0);
    EXPECT_DOUBLE_EQ(0.0, solver.maxDensityErrorRatio());

    solver.setMaxNumberOfIterations(10);
    EXPECT_DOUBLE_EQ(10, solver.maxNumberOfIterations());
}
