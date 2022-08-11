//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.force/binary_file_reader_writer.h"

namespace vox::flex {
/** \brief Class to manage the current simulation time and the time step size.
 * This class is a singleton.
 */
class TimeManager {
private:
    double time;
    static TimeManager *current;
    double h;

public:
    TimeManager();
    ~TimeManager();

    // Singleton
    static TimeManager *getCurrent();
    static void setCurrent(TimeManager *tm);
    static bool hasCurrent();

    [[nodiscard]] double getTime() const;
    void setTime(double t);
    [[nodiscard]] double getTimeStepSize() const;
    void setTimeStepSize(double tss);

    void saveState(BinaryFileWriter &binWriter) const;
    void loadState(BinaryFileReader &binReader);
};
}  // namespace vox::flex
