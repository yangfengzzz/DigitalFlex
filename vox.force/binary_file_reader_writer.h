//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <fstream>
#include <iostream>

#include "vox.math/matrix.h"

namespace vox::flex {
//MARK: - BinaryFileWriter
class BinaryFileWriter {
public:
    std::ofstream m_file;

public:
    bool openFile(const std::string &fileName) {
        m_file.open(fileName, std::ios::out | std::ios::binary);
        if (!m_file.is_open()) {
            std::cout << "Cannot open file.\n";
            return false;
        }
        return true;
    }

    void closeFile() { m_file.close(); }

    void writeBuffer(const char *buffer, size_t size) { m_file.write(buffer, size); }

    template <typename T>
    void write(const T &v) {
        writeBuffer((char *)&v, sizeof(T));
    }

    void write(const std::string &str) {
        write((unsigned int)str.size());
        writeBuffer(str.c_str(), str.size());
    }

    template <typename T>
    void writeMatrix(const T &m) {
        writeBuffer((char *)m.data(), m.size() * sizeof(m.data()[0]));
    }

    template <typename T, int Rows, int Cols>
    void writeMatrixX(const Matrix<T, Rows, Cols> &m) {
        const size_t rows = m.rows();
        const size_t cols = m.cols();
        write(rows);
        write(cols);

        writeBuffer((char *)m.data(), rows * cols * sizeof(T));
    }

    template <typename T>
    void writeVector(const std::vector<T> &m) {
        write(m.size());
        writeBuffer((char *)m.data(), m.size() * sizeof(T));
    }
};

//MARK: - BinaryFileReader
class BinaryFileReader {
public:
    std::ifstream m_file;

public:
    bool openFile(const std::string &fileName) {
        m_file.open(fileName, std::ios::in | std::ios::binary);
        if (!m_file.is_open()) {
            std::cout << "Cannot open file.\n";
            return false;
        }
        return true;
    }

    void closeFile() { m_file.close(); }

    void readBuffer(char *buffer, size_t size) { m_file.read(buffer, size); }

    template <typename T>
    void read(T &v) {
        readBuffer((char *)&v, sizeof(T));
    }

    void read(std::string &str) {
        unsigned int len;
        read(len);
        char *temp = new char[len + 1u];
        readBuffer(temp, len);
        temp[len] = '\0';
        str = std::string(temp);
        delete[] temp;
    }

    template <typename T>
    void readMatrix(T &m) {
        readBuffer((char *)m.data(), m.size() * sizeof(m.data()[0]));
    }

    template <typename T, int Rows, int Cols>
    void readMatrixX(Matrix<T, Rows, Cols> &m) {
        size_t rows, cols;
        read(rows);
        read(cols);
        m.resize(rows, cols);

        readBuffer((char *)m.data(), rows * cols * sizeof(T));
    }

    template <typename T>
    void readVector(std::vector<T> &m) {
        size_t size;
        read(size);
        m.resize(size);
        readBuffer((char *)m.data(), m.size() * sizeof(T));
    }
};

}  // namespace vox::flex
