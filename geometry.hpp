#ifndef DYNEARTHSOL3D_GEOMETRY_HPP
#define DYNEARTHSOL3D_GEOMETRY_HPP

void compute_volume(const double2d &coord, const int2d &connectivity,
                    double_vec &volume, double_vec &volume_n);

double compute_dt(const Param& param, const Variables& var);

void compute_mass(const Param &param,
                  const double2d &coord, const int2d &connectivity,
                  const double_vec &volume, const MatProps &mat,
                  double_vec &mass, double_vec &tmass);

void compute_shape_fn(const double2d &coord, const int2d &connectivity,
                      const double_vec &volume,
                      const std::vector<int_vec> &egroups,
                      double2d &shpdx, double2d &shpdy, double2d &shpdz);

#endif
