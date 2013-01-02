#ifndef DYNEARTHSOL3D_RHEOLOGY_HPP
#define DYNEARTHSOL3D_RHEOLOGY_HPP

void update_stress(const Variables& var, tensor_t& stress,
                   tensor_t& strain, double_vec& plstrain,
                   tensor_t& strain_rate);

#endif
