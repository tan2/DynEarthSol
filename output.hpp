#ifndef DYNEARTHSOL3D_OUTPUT_HPP
#define DYNEARTHSOL3D_OUTPUT_HPP

#include "array2d.hpp"

static const char revision_str[] = "# DynEarthSol ndims="
#ifdef THREED
    "3"
#else
    "2"
#endif
    " revision=1\n";

class Output
{
private:
    const std::string &modelname;
    const double start_time;
    const int average_interval;
    int frame;

    // stuffs for averging fields
    double time0;
    array_t coord0;
    tensor_t strain0;
    tensor_t stress_avg;
    double_vec delta_plstrain_avg;

    // stuffs for header
    long eof_pos;
    char *header;
    char *hd_pos;
    static const int headerlen = 4096;

    template <typename T>
    void write_array(std::FILE* f, const std::vector<T>& A, const char *name);

    template <typename T, int N>
    void write_array(std::FILE* f, const Array2D<T,N>& A, const char *name);

    void write_header(std::FILE* f, const char *name);
    void write_info(const Variables& var, double dt);

public:
    Output(const Param& param, double start_time, int start_frame);
    ~Output();
    void write(const Variables& var, bool is_averaged=true);
    void average_fields(Variables& var);

};


void restart();

#endif
