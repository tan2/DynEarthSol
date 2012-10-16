#ifndef DYNEARTHSOL3D_PARAMETERS_HPP
#define DYNEARTHSOL3D_PARAMETERS_HPP

#include <string>
#include <boost/multi_array.hpp>

typedef boost::multi_array<double,2> double2d;
typedef boost::multi_array<int,2> int2d;
typedef boost::multi_array<int,1> int1d;

typedef boost::multi_array_ref<double,2> double2d_ref;
typedef boost::multi_array_ref<int,2> int2d_ref;
typedef boost::multi_array_ref<int,1> int1d_ref;

//
// Structures for input parameters
//
struct Sim {
    double max_time;
    int max_steps;
    int output_step_interval;
    double output_time_interval;
    bool is_restarting;

    std::string modelname;
};

struct Mesh {
    double xlength, ylength, zlength;
    double resolution;
};

struct Param {
    Sim sim;
    Mesh mesh;
};


//
// Structures for model variables
//
struct Variables {
    double time;
    double dt;
    int steps;
    int frame;

    int nnode;
    int nelem;
    int nseg;


    // These 4 arrays are allocated by external library
    double2d_ref* coord;
    int2d_ref* connectivity;
    int2d_ref* segment;
    int1d_ref* segflag;
};

#endif
