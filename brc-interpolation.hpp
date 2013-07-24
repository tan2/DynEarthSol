#ifndef DYNEARTHSOL3D_BRC_INTERPOLATION_HPP
#define DYNEARTHSOL3D_BRC_INTERPOLATION_HPP

void barycentric_node_interpolation(Variables &var,
                                    const array_t &old_coord,
                                    const conn_t &old_connectivity);

#endif
