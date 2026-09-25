function solve_bleb_population(mesh_dir, specs)
% SOLVE_BLEB_POPULATION  For every bleb_XXX_mesh.mat in mesh_dir: rebuild the boundary,
% solve the BEM problem, and APPEND the compact solution to the same file as `bem_sol`
% (e, h, k0 + nelem/pos_checksum for the mismatch check on reload). Files that already
% contain bem_sol are skipped, so reruns resume.
%
%   specs = make_imaging_specs();
%   solve_bleb_population(mesh_dir, specs);

    if ~isfield(specs, 'tau_method'), specs.tau_method = 'native'; end   % tolerate a stale specs struct
    sys   = build_imaging_system(specs);
    einc  = sys.einc;
    files = dir(fullfile(mesh_dir, 'bleb_*_mesh.mat'));
    nf    = numel(files);
    fprintf('Solving %d blebs...\n', nf);

    for f = 11:50
        fn = fullfile(mesh_dir, files(f).name);
        fprintf('Solving bleb number %d ...\n', f);
        if ismember('bem_sol', who('-file', fn))
            continue                                   % already solved
        end

        M   = load(fn, 'mesh');
        tau = rebuild_bleb_tau(M.mesh, sys.mat_set, specs.idx_water, specs.idx_lipid, ...
                               specs.idx_aqueous, specs.tau_method);

        qinc = einc(tau, 'layer', sys.layer);          % incident field incl. substrate reflection
        bem  = stratified.bemsolver(tau, sys.layer, 'order', specs.bem_order);
        sol  = bem \ qinc;

        bem_sol = struct('e', sol.e, 'h', sol.h, 'k0', sol.k0, ...
                         'nelem', numel(tau), ...
                         'pos_checksum', sum(vertcat(tau.pos), 'all'), ...
                         'lambda', specs.lambda, 'pol', specs.pol, 'dir', specs.dir);
        save(fn, 'bem_sol', '-append');
        fprintf('  [%d/%d] solved %s (%d elements)\n', f, nf, files(f).name, numel(tau));
    end
end
