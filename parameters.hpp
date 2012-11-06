#ifndef DYNEARTHSOL3D_PARAMETERS_HPP
#define DYNEARTHSOL3D_PARAMETERS_HPP

#include <string>
#include <utility>
#include <vector>
#include <boost/multi_array.hpp>

typedef std::pair<double,double> double_pair;

typedef std::vector<double> double_vec;
typedef std::vector<int> int_vec;

typedef boost::multi_array<double,2> double2d;
typedef boost::multi_array<double,1> double1d;
typedef boost::multi_array<int,2> int2d;
typedef boost::multi_array<int,1> int1d;

typedef boost::multi_array_ref<double,2> double2d_ref;
typedef boost::multi_array_ref<int,2> int2d_ref;
typedef boost::multi_array_ref<int,1> int1d_ref;

//
// Structures for input parameters
//
struct Sim {
    double max_time_in_yr;
    double output_time_interval_in_yr;
    int max_steps;
    int output_step_interval;
    bool is_restarting;

    std::string modelname;
};

struct Mesh {
    int meshing_option;

    double xlength, ylength, zlength;
    double resolution;
    // for 2D only
    double min_angle;
    // for 3D only
    double min_tet_angle, max_ratio;

    double_pair refined_zonex, refined_zoney, refined_zonez;
};

struct Param {
    Sim sim;
    Mesh mesh;

    double surface_temperature;
    double mantle_temperature;

    double strain_inert;
    double maxvbcval;
};


//
// Structures for model variables
//
class MatProps;
struct Variables {
    double time;
    double dt;
    int steps;
    int frame;

    int nnode;
    int nelem;
    int nseg;

    // These 4 arrays are allocated by external library
    double2d_ref *coord;
    int2d_ref *connectivity;
    int2d_ref *segment;
    int1d_ref *segflag;

    int_vec *bcflag;

    std::vector<int_vec> *support;

    double_vec *volume, *volume_old, *volume_n;
    double_vec *mass, *tmass;
    double_vec *jacobian, *ejacobian;
    double_vec *temperature, *plstrain;
    double_vec *tmp0;

    double2d *vel, *force;
    double2d *strain_rate, *strain, *stress;
    double2d *shpdx, *shpdy, *shpdz;

    MatProps *mat;
};

#endif
