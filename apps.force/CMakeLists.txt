cmake_minimum_required(VERSION 3.12)

project(apps.force LANGUAGES C CXX)

file(GLOB sources
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(${PROJECT_NAME} ${sources})

target_include_directories(${PROJECT_NAME} PUBLIC ../
        ${CMAKE_SOURCE_DIR}/third_party/googletest/googlemock/include
        ${CMAKE_SOURCE_DIR}/third_party/googletest/googletest/include
        ${CMAKE_SOURCE_DIR}/third_party/pystring
        ${CMAKE_SOURCE_DIR}/third_party/cnpy)

# Link third party libraries
target_link_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/third_party/googletest/build/lib
        ${CMAKE_SOURCE_DIR}/third_party/cnpy/build
        ${CMAKE_SOURCE_DIR}/third_party/pystring/build)
target_link_libraries(${PROJECT_NAME} PUBLIC spdlog vox.base vox.math vox.geometry vox.force
        libgmock.a libgmock_main.a libgtest.a libgtest_main.a cnpy pystring)