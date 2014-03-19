#ifndef DYNEARTHSOL3D_BINARYIO_HPP
#define DYNEARTHSOL3D_BINARYIO_HPP

#include "array2d.hpp"


class BinaryOutput
{
private:
    long eof_pos;
    char *header;
    char *hd_pos;
    std::FILE* f;

    void write_header(const char *name);

public:
    BinaryOutput(const char *filename);
    ~BinaryOutput();

    template <typename T>
    void write_array(const std::vector<T>& A, const char *name);

    template <typename T, int N>
    void write_array(const Array2D<T,N>& A, const char *name);
};

#endif
