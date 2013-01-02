#ifndef DYNEARTHSOL3D_GEOMETRY_HPP
#define DYNEARTHSOL3D_GEOMETRY_HPP

void compute_volume(const array_t &coord, const conn_t &connectivity,
                    const std::vector<int_vec> &egroups,
                    double_vec &volume, double_vec &volume_n);

double compute_dt(const Param& param, const Variables& var);

void compute_mass(const Param &param,
                  const std::vector<int_vec> &egroups, const conn_t &connectivity,
                  const double_vec &volume, const MatProps &mat,
                  double max_vbc_val,
                  double_vec &mass, double_vec &tmass);

void compute_shape_fn(const array_t &coord, const conn_t &connectivity,
                      const double_vec &volume,
                      const std::vector<int_vec> &egroups,
                      shapefn &shpdx, shapefn &shpdy, shapefn &shpdz);

#endif
