function build_sphere_population(mesh_dir, N, d_mean, d_sd, specs, seed, prefix)
% BUILD_SPHERE_POPULATION  Generate N homogeneous-sphere meshes (full or empty LNPs) and
% save each as its own lean file sphere_XXX_mesh.mat in mesh_dir, same pattern as the
% blebbed-LNP population builder: one small file per particle, resumable, seeded.
%
%   build_sphere_population(mesh_dir, N, d_mean, d_sd, specs, seed, 'full')
%   build_sphere_population(mesh_dir, N, d_mean, d_sd, specs, seed, 'empty')
%
%   mesh_dir       : output folder (created if missing)
%   N              : number of particles
%   d_mean, d_sd   : diameter distribution (nm), normal, resampled if a draw is <= 0
%   specs          : make_imaging_specs() output (uses specs.nverts, specs.gap)
%   seed           : RNG seed, so the population is exactly reproducible
%   prefix         : 'full' or 'empty' -- just a label used in the saved params and in the
%                    filename prefix ('full_001_mesh.mat', 'empty_001_mesh.mat', ...), so
%                    the two populations can live in the same or different folders without
%                    colliding.
%
% Each saved mesh.params also carries d_mean/d_sd/seed/index/label, so any particle's file
% alone tells you exactly how it was generated.

    if ~exist(mesh_dir, 'dir'), mkdir(mesh_dir); end
    rng(seed);

    for k = 1:N
        fn = fullfile(mesh_dir, sprintf('%s_%03d_mesh.mat', prefix, k));

        d = -1;
        while d <= 0                 % reject non-physical draws from the normal tail
            d = d_mean + d_sd * randn;
        end

        if exist(fn, 'file'), continue; end   % resume: RNG stream still advances above

        mesh = build_sphere_mesh(d, specs.nverts, specs.gap);
        mesh.params.label  = prefix;
        mesh.params.d_mean = d_mean;
        mesh.params.d_sd   = d_sd;
        mesh.params.seed   = seed;
        mesh.params.index  = k;

        save(fn, 'mesh');
        fprintf('saved %s (d = %.1f nm)\n', fn, d);
    end
end
