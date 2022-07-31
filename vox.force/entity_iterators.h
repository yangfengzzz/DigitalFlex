//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <array>
#include <iterator>

#include "vox.force/half_edge.h"
#include "vox.math/vector3.h"

namespace vox::flex {

class TriangleMesh;

class FaceContainer;
class FaceIterator {
public:
    typedef std::array<unsigned int, 3> value_type;
    typedef ptrdiff_t difference_type;
    typedef value_type *pointer;
    typedef value_type &reference;
    typedef std::random_access_iterator_tag iterator_category;
    typedef FaceIterator Self;

    FaceIterator() = delete;

    reference operator*();

    bool operator<(Self const &other) const { return m_index < other.m_index; }
    bool operator==(Self const &other) const { return m_index == other.m_index; }

    bool operator!=(Self const &other) const { return !(*this == other); }

    inline Self &operator++() {
        ++m_index;
        return *this;
    }
    inline Self &operator--() {
        --m_index;
        return *this;
    }

    inline Self operator+(Self const &rhs) { return {m_index + rhs.m_index, m_mesh}; }
    inline difference_type operator-(Self const &rhs) const { return m_index - rhs.m_index; }
    inline Self operator-(int const &rhs) { return {m_index - rhs, m_mesh}; }

    [[nodiscard]] unsigned int vertex(unsigned int i) const;
    unsigned int &vertex(unsigned int i);

private:
    friend class FaceContainer;
    FaceIterator(unsigned int index, TriangleMesh *mesh) : m_index(index), m_mesh(mesh) {}

    unsigned int m_index;
    TriangleMesh *m_mesh;
};

class FaceConstIterator {
public:
    typedef std::array<unsigned int, 3> const value_type;
    typedef ptrdiff_t difference_type;
    typedef value_type *pointer;
    typedef value_type &reference;
    typedef std::random_access_iterator_tag iterator_category;
    typedef FaceConstIterator Self;

    FaceConstIterator() = delete;

    reference operator*();

    bool operator<(Self const &other) const { return m_index < other.m_index; }
    bool operator==(Self const &other) const { return m_index == other.m_index; }

    bool operator!=(Self const &other) const { return !(*this == other); }

    inline Self &operator++() {
        ++m_index;
        return *this;
    }
    inline Self &operator--() {
        --m_index;
        return *this;
    }

    inline Self operator+(Self const &rhs) const { return {m_index + rhs.m_index, m_mesh}; }
    inline difference_type operator-(Self const &rhs) const { return m_index - rhs.m_index; }
    inline Self operator-(int const &rhs) const { return {m_index - rhs, m_mesh}; }

    [[nodiscard]] unsigned int vertex(unsigned int i) const;
    unsigned int &vertex(unsigned int i);

private:
    friend class FaceConstContainer;
    FaceConstIterator(unsigned int index, TriangleMesh const *mesh) : m_index(index), m_mesh(mesh) {}

    unsigned int m_index;
    TriangleMesh const *m_mesh;
};

class IncidentFaceContainer;
class IncidentFaceIterator {
public:
    typedef Halfedge value_type;
    typedef ptrdiff_t difference_type;
    typedef value_type *pointer;
    typedef value_type &reference;
    typedef std::random_access_iterator_tag iterator_category;
    typedef IncidentFaceIterator Self;

    value_type operator*() { return m_h; }
    Self &operator++();
    bool operator==(Self const &other) const { return m_h == other.m_h; }

    bool operator!=(Self const &other) const { return !(*this == other); }

private:
    friend class IncidentFaceContainer;
    IncidentFaceIterator(unsigned int v, TriangleMesh const *mesh);
    IncidentFaceIterator() : m_h(), m_begin(), m_mesh(nullptr) {}

    Halfedge m_h, m_begin;
    TriangleMesh const *m_mesh;
};

class VertexContainer;
class VertexIterator {
public:
    typedef Vector3D value_type;
    typedef ptrdiff_t difference_type;
    typedef value_type *pointer;
    typedef value_type &reference;
    typedef std::random_access_iterator_tag iterator_category;
    typedef VertexIterator Self;

    VertexIterator() = delete;

    reference operator*();

    bool operator<(Self const &other) const { return m_index < other.m_index; }
    bool operator==(Self const &other) const { return m_index == other.m_index; }

    bool operator!=(Self const &other) const { return !(*this == other); }

    inline Self &operator++() {
        ++m_index;
        return *this;
    }
    inline Self &operator--() {
        --m_index;
        return *this;
    }

    inline Self operator+(Self const &rhs) const { return {m_index + rhs.m_index, m_mesh}; }
    inline difference_type operator-(Self const &rhs) const { return m_index - rhs.m_index; }
    inline Self operator-(int const &rhs) const { return {m_index - rhs, m_mesh}; }

    [[nodiscard]] unsigned int index() const;

private:
    friend class VertexContainer;
    VertexIterator(unsigned int index, TriangleMesh *mesh) : m_index(index), m_mesh(mesh) {}

    unsigned int m_index;
    TriangleMesh *m_mesh;
};

class VertexConstContainer;
class VertexConstIterator {
public:
    typedef Vector3D const value_type;
    typedef ptrdiff_t difference_type;
    typedef value_type *pointer;
    typedef value_type &reference;
    typedef std::random_access_iterator_tag iterator_category;
    typedef VertexConstIterator Self;

    VertexConstIterator() = delete;

    reference operator*();

    bool operator<(Self const &other) const { return m_index < other.m_index; }
    bool operator==(Self const &other) const { return m_index == other.m_index; }

    bool operator!=(Self const &other) const { return !(*this == other); }

    inline Self &operator++() {
        ++m_index;
        return *this;
    }
    inline Self &operator--() {
        --m_index;
        return *this;
    }

    inline Self operator+(Self const &rhs) const { return {m_index + rhs.m_index, m_mesh}; }
    inline difference_type operator-(Self const &rhs) const { return m_index - rhs.m_index; }
    inline Self operator-(int const &rhs) const { return {m_index - rhs, m_mesh}; }

    [[nodiscard]] unsigned int index() const;

private:
    friend class VertexConstContainer;
    VertexConstIterator(unsigned int index, TriangleMesh const *mesh) : m_index(index), m_mesh(mesh) {}

    unsigned int m_index;
    TriangleMesh const *m_mesh;
};
}  // namespace vox::flex
