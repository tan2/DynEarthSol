#ifndef DYNEARTHSOL3D_MESH_HPP
#define DYNEARTHSOL3D_MESH_HPP

void points_to_mesh(const Param &param, Variables &var,
                    int npoints, double *points,
                    int n_init_segments, int *init_segments, int *init_segflags,
                    double max_elem_size);

void create_new_mesh(const Param&, Variables&);

#endif
