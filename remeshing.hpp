#ifndef DYNEARTHSOL3D_REMESHING_HPP
#define DYNEARTHSOL3D_REMESHING_HPP

bool bad_mesh_quality(const Param&, const Variables&, int&);
void remesh(const Param&, Variables&);

#endif
