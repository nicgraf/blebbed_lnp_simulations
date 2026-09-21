function [particle, sol] = load_bleb_solution(fn, sys)
% LOAD_BLEB_SOLUTION  Load one bleb_XXX_mesh.mat (with bem_sol appended by
% solve_bleb_population), rebuild tau, verify it matches the solved mesh, and return
%   particle : struct with .params (d_core, d_bleb, neck_overlap, theta/phi, ...) and .tau
%   sol      : stratified.solution ready for propagate_iscat_image / _zstack
% `sys` is the output of build_imaging_system (or specs, which is then built here).

    if ~isfield(sys, 'tau_method'), sys.tau_method = 'native'; end       % tolerate a stale specs/sys struct
    if ~isfield(sys, 'lens'), sys = build_imaging_system(sys); end

    M = load(fn, 'mesh', 'bem_sol');
    assert(isfield(M, 'bem_sol'), 'No bem_sol in %s -- run solve_bleb_population first.', fn);

    tau = rebuild_bleb_tau(M.mesh, sys.mat_set, sys.idx_water, sys.idx_lipid, ...
                           sys.idx_aqueous, sys.tau_method);
    B = M.bem_sol;
    assert(numel(tau) == B.nelem && ...
           abs(sum(vertcat(tau.pos), 'all') - B.pos_checksum) <= 1e-9 * max(abs(B.pos_checksum), 1), ...
           'Rebuilt mesh does not match the saved solution in %s.', fn);

    sol = stratified.solution(tau, B.k0, B.e, B.h);
    sol.layer = sys.layer;

    particle = struct('params', M.mesh.params, 'tau', tau);
end
