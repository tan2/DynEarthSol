#ifndef DYNEARTHSOL3D_RHEOLOGY_HPP
#define DYNEARTHSOL3D_RHEOLOGY_HPP

void update_stress(const Variables& var, double2d& stress,
                   double2d& strain, double_vec& plstrain,
                   double2d& strain_rate);

#endif
