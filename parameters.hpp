#ifndef DYNEARTHSOL3D_PARAMETERS_HPP
#define DYNEARTHSOL3D_PARAMETERS_HPP

struct Sim {
    double max_time;
    int max_steps;
    int output_step_interval;
    double output_time_interval;
    bool is_restarting;
};

struct Mesh {
    double xlength, ylength, zlength;
};

struct Param {
    Sim sim;
    Mesh mesh;
};

#endif
