//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "time_manager.h"

namespace vox::flex {

TimeManager* TimeManager::current = nullptr;

TimeManager::TimeManager() {
    time = 0;
    h = 0.001;
}

TimeManager::~TimeManager() { current = nullptr; }

TimeManager* TimeManager::getCurrent() {
    if (current == nullptr) {
        current = new TimeManager();
    }
    return current;
}

void TimeManager::setCurrent(TimeManager* tm) { current = tm; }

bool TimeManager::hasCurrent() { return (current != nullptr); }

double TimeManager::getTime() const { return time; }

void TimeManager::setTime(double t) { time = t; }

double TimeManager::getTimeStepSize() const { return h; }

void TimeManager::setTimeStepSize(double tss) { h = tss; }

void TimeManager::saveState(BinaryFileWriter& binWriter) const {
    binWriter.write(time);
    binWriter.write(h);
}

void TimeManager::loadState(BinaryFileReader& binReader) {
    binReader.read(time);
    binReader.read(h);
}
}  // namespace vox::flex