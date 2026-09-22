function [particle_out, sol] = load_sphere_solution(fn, sys, prefix)
% LOAD_SPHERE_SOLUTION  Load one <prefix>_XXX_mesh.mat (with bem_sol appended by
% solve_sphere_population), rebuild tau, verify it matches the solved mesh, and return
%   particle_out : struct with .params (d, nverts, gap, label, ...) and .tau
%   sol          : stratified.solution ready for propagate_iscat_image / _zstack
%
%   [particle_out, sol] = load_sphere_solution(fn, sys, 'full')
%
% `sys` is the output of build_imaging_system (or specs, which is then built here).

    if ~isfield(sys, 'lens'), sys = build_imaging_system(sys); end

    switch prefix
        case 'full',  idx_particle = sys.idx_full;
        case 'empty', idx_particle = sys.idx_empty;
        otherwise
            error('load_sphere_solution:prefix', 'prefix must be ''full'' or ''empty''.');
    end

    M = load(fn, 'mesh', 'bem_sol');
    assert(isfield(M, 'bem_sol'), 'No bem_sol in %s -- run solve_sphere_population first.', fn);

    tau = rebuild_sphere_tau(M.mesh, sys.mat_set, idx_particle, sys.idx_water);
    B = M.bem_sol;
    assert(numel(tau) == B.nelem && ...
           abs(sum(vertcat(tau.pos), 'all') - B.pos_checksum) <= 1e-9 * max(abs(B.pos_checksum), 1), ...
           'Rebuilt mesh does not match the saved solution in %s.', fn);

    sol = stratified.solution(tau, B.k0, B.e, B.h);
    sol.layer = sys.layer;

    particle_out = struct('params', M.mesh.params, 'tau', tau);
end
