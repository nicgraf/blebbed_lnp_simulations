function solve_sphere_population(mesh_dir, specs, prefix, varargin)
% SOLVE_SPHERE_POPULATION  For every <prefix>_XXX_mesh.mat in mesh_dir: rebuild the
% boundary, solve the BEM problem, and APPEND the compact solution to the same file as
% `bem_sol` (e, h, k0 + nelem/pos_checksum for the mismatch check on reload). Files that
% already contain bem_sol are skipped, so reruns resume. Same storage/solve pattern as
% solve_bleb_population.m, generalized to the single-boundary sphere case.
%
%   specs = make_imaging_specs();
%   solve_sphere_population(mesh_dir, specs, 'full');                        % parallel
%   solve_sphere_population(mesh_dir, specs, 'empty', 'UseParallel', false); % serial
%
% prefix : 'full' or 'empty' -- selects specs.idx_full/idx_empty and the mat_set entry to
% solve against, and which <prefix>_XXX_mesh.mat files to pick up.

    p = inputParser;
    addParameter(p, 'UseParallel', true, @islogical);
    addParameter(p, 'NumWorkers', []);
    addParameter(p, 'bem_order', []);   % [] = use specs.bem_order; spheres often converge
                                        % at a lower order than the fused bleb mesh
    parse(p, varargin{:});
    useParallel = p.Results.UseParallel && license('test', 'Distrib_Computing_Toolbox') ...
                  && ~isempty(which('parpool'));
    if p.Results.UseParallel && ~useParallel
        warning('solve_sphere_population:noPCT', ...
            'Parallel Computing Toolbox not available -- running serially.');
    end

    switch prefix
        case 'full',  idx_particle = specs.idx_full;
        case 'empty', idx_particle = specs.idx_empty;
        otherwise
            error('solve_sphere_population:prefix', 'prefix must be ''full'' or ''empty''.');
    end
    bemOrder = p.Results.bem_order;
    if isempty(bemOrder), bemOrder = specs.bem_order; end

    sys   = build_imaging_system(specs);
    einc  = sys.einc;
    files = dir(fullfile(mesh_dir, [prefix '_*_mesh.mat']));
    nf    = numel(files);
    fprintf('Solving %d %s spheres (%s)...\n', nf, prefix, ternary(useParallel, 'parallel', 'serial'));

    if useParallel
        pool = gcp('nocreate');
        if isempty(pool)
            try
                % see solve_bleb_population.m for why JobStorageLocation is forced here
                c = parcluster('local');
                c.JobStorageLocation = tempdir;
                if ~isempty(p.Results.NumWorkers)
                    pool = parpool(c, p.Results.NumWorkers);
                else
                    pool = parpool(c);
                end
            catch ME
                warning('solve_sphere_population:parpoolFailed', ...
                    'parpool failed to start (%s) -- falling back to serial.', ME.message);
                useParallel = false;
            end
        end
    end

    if useParallel
        clientPath = path;
        pctRunOnAll(sprintf('addpath(''%s'');', strrep(clientPath, '''', '''''')));
        parfor f = 1:nf
            solve_one(mesh_dir, files(f).name, f, nf, sys, einc, idx_particle, specs, bemOrder); %#ok<PFBNS>
        end
    else
        for f = 1:nf
            solve_one(mesh_dir, files(f).name, f, nf, sys, einc, idx_particle, specs, bemOrder);
        end
    end
end

function solve_one(mesh_dir, fname, f, nf, sys, einc, idx_particle, specs, bemOrder)
    fn = fullfile(mesh_dir, fname);
    if ismember('bem_sol', who('-file', fn))
        fprintf('  [%d/%d] %s already solved, skipping\n', f, nf, fname);
        return
    end

    M   = load(fn, 'mesh');
    tau = rebuild_sphere_tau(M.mesh, sys.mat_set, idx_particle, sys.idx_water);

    qinc = einc(tau, 'layer', sys.layer);
    bem  = stratified.bemsolver(tau, sys.layer, 'order', bemOrder);
    sol  = bem \ qinc;

    bem_sol = struct('e', sol.e, 'h', sol.h, 'k0', sol.k0, ...
                     'nelem', numel(tau), ...
                     'pos_checksum', sum(vertcat(tau.pos), 'all'), ...
                     'lambda', specs.lambda, 'pol', specs.pol, 'dir', specs.dir);
    save(fn, 'bem_sol', '-append');
    fprintf('  [%d/%d] solved %s (%d elements)\n', f, nf, fname, numel(tau));
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
