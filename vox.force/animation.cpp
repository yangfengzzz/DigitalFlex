// Copyright (c) 2022 Feng Yang
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.force/animation.h"

#include "vox.base/logging.h"
#include "vox.base/timer.h"

using namespace vox;

Animation::Animation() = default;

Animation::~Animation() = default;

void Animation::update(const Frame& frame) {
    utility::Timer timer;

    LOGI("Begin updating frame: {} timeIntervalInSeconds: {}  (1/{}) seconds", frame.index, frame.timeIntervalInSeconds,
         1.0 / frame.timeIntervalInSeconds)

    onUpdate(frame);

    LOGI("End updating frame (took {} seconds)", timer.Elapsed())
}
