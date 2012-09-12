#ifndef DYNEARTHSOL3D_PARAMETERS_HPP
#define DYNEARTHSOL3D_PARAMETERS_HPP

struct Sim {
    int max_steps;
    double max_time;
};

struct Mesh {
    double xlength, ylength, zlength;
};

struct Param {
    Sim sim;
    Mesh mesh;
};

#endif
