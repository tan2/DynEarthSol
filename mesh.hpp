#ifndef DYNEARTHSOL3D_MESH_HPP
#define DYNEARTHSOL3D_MESH_HPP

void points_to_new_mesh(const Mesh &mesh, int npoints, const double *points,
                        int n_init_segments, const int *init_segments, const int *init_segflags,
                        int n_regions, const double *regattr,
                        double max_elem_size, int vertex_per_polygon,
                        int &nnode, int &nelem, int &nseg,
                        double *&pcoord, int *&pconnectivity,
                        int *&psegment, int *&psegflag, double *&pregattr);
void points_to_mesh(const Param &param, Variables &var,
                    int npoints, const double *points,
                    int n_init_segments, const int *init_segments, const int *init_segflags, const double *regattr,
                    double max_elem_size, int vertex_per_polygon);
void create_boundary_flags2(uint_vec &bcflag, const segment_t &segment,
                            const segflag_t &segflag);
void create_boundary_flags(Variables& var);
void create_boundary_nodes(Variables& var);
void create_boundary_facets(Variables& var);
void create_support(Variables& var);
void create_elem_groups(Variables& var);
void create_elemmarkers(const Param&, Variables&);
void create_markers(const Param&, Variables&);
void create_new_mesh(const Param&, Variables&);

#endif
