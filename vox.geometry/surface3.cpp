// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include "vox.geometry/surface3.h"

#include <algorithm>

using namespace vox;

Surface3::Surface3(const Transform3D &transform_, bool isNormalFlipped_)
    : transform(transform_), isNormalFlipped(isNormalFlipped_) {}

Surface3::Surface3(const Surface3 &other) : transform(other.transform), isNormalFlipped(other.isNormalFlipped) {}

Surface3::~Surface3() {}

Point3D Surface3::closestPoint(const Point3D &otherPoint) const {
    return transform.toWorld(closestPointLocal(transform.toLocal(otherPoint)));
}

BoundingBox3D Surface3::boundingBox() const { return transform.toWorld(boundingBoxLocal()); }

bool Surface3::intersects(const Ray3D &ray) const { return intersectsLocal(transform.toLocal(ray)); }

double Surface3::closestDistance(const Point3D &otherPoint) const {
    return closestDistanceLocal(transform.toLocal(otherPoint));
}

SurfaceRayIntersection3 Surface3::closestIntersection(const Ray3D &ray) const {
    auto result = closestIntersectionLocal(transform.toLocal(ray));
    result.point = transform.toWorld(result.point);
    result.normal = transform.toWorldDirection(result.normal);
    result.normal *= (isNormalFlipped) ? -1.0 : 1.0;
    return result;
}

Vector3D Surface3::closestNormal(const Point3D &otherPoint) const {
    auto result = transform.toWorldDirection(closestNormalLocal(transform.toLocal(otherPoint)));
    result *= (isNormalFlipped) ? -1.0 : 1.0;
    return result;
}

bool Surface3::intersectsLocal(const Ray3D &rayLocal) const {
    auto result = closestIntersectionLocal(rayLocal);
    return result.isIntersecting;
}

void Surface3::updateQueryEngine() {
    // Do nothing
}

bool Surface3::isBounded() const { return true; }

bool Surface3::isValidGeometry() const { return true; }

bool Surface3::isInside(const Point3D &otherPoint) const {
    return isNormalFlipped == !isInsideLocal(transform.toLocal(otherPoint));
}

double Surface3::closestDistanceLocal(const Point3D &otherPointLocal) const {
    return otherPointLocal.distanceTo(closestPointLocal(otherPointLocal));
}

bool Surface3::isInsideLocal(const Point3D &otherPointLocal) const {
    Point3D cpLocal = closestPointLocal(otherPointLocal);
    Vector3D normalLocal = closestNormalLocal(otherPointLocal);
    return (otherPointLocal - cpLocal).dot(normalLocal) < 0.0;
}
