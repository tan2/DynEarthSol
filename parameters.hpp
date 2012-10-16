#ifndef DYNEARTHSOL3D_PARAMETERS_HPP
#define DYNEARTHSOL3D_PARAMETERS_HPP

#include <string>
#include <boost/multi_array.hpp>

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
    boost::multi_array_ref<double,2>* coord;
    boost::multi_array_ref<int,2>* connectivity;
    boost::multi_array_ref<int,2>* segment;
    boost::multi_array_ref<int,1>* segflag;
};

#endif
