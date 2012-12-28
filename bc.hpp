#ifndef DYNEARTHSOL3D_BC_HPP
#define DYNEARTHSOL3D_BC_HPP

bool is_on_boundary(const Variables &var, int node);
double find_max_vbc(const BC &bc);
void apply_vbcs(const Param &param, const Variables &var, arrayd2 &vel);
void apply_stress_bcs(const Param& param, const Variables& var, arrayd2& force);

#endif
