function tau = rebuild_sphere_tau(mesh, mat_set, idx_particle, idx_water)
% REBUILD_SPHERE_TAU  Turn a lean mesh struct (from build_sphere_mesh) back into the
% single-boundary `tau` that the BEM solver / stratified.solution need.
%
%   tau = rebuild_sphere_tau(mesh, mat_set, idx_particle, idx_water)
%
% idx_particle is idx_full or idx_empty (from make_imaging_specs), i.e. which material in
% mat_set this sphere is made of; idx_water is the surrounding medium.

    p = particle(mesh.node, double(mesh.faces));
    tau = BoundaryEdge(mat_set, p, [idx_particle, idx_water]);   % particle inside, water outside
end
