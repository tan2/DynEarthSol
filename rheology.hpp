#ifndef DYNEARTHSOL3D_RHEOLOGY_HPP
#define DYNEARTHSOL3D_RHEOLOGY_HPP

void update_stress(const Variables& var, tensord2& stress,
                   tensord2& strain, double_vec& plstrain,
                   tensord2& strain_rate);

#endif
