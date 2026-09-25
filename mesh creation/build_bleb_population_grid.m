function build_bleb_population_grid(mesh_dir, specs, varargin)
% BUILD_BLEB_POPULATION_GRID  Deterministic (d_core, d_bleb) size-sweep of blebbed-LNP
% meshes, for rigorous/systematic testing -- every combination on a regular grid,
% not a randomly-sampled population. One mesh per grid cell, saved as its own lean
% file, same one-file-per-particle / resumable pattern as build_sphere_population.m
% and the original blebbed-LNP loop in lnp_mesh_creation.mlx.
%
%   build_bleb_population_grid(mesh_dir, specs)
%
%       -> 22 core sizes x 22 bleb sizes, each 10-120 nm (linspace), fixed
%          OverlapFraction = 0.6, fixed ThetaTiltDeg = 90 / PhiAzimuthDeg = 0
%          ("flat on the glass" -- see below) = 484 particles.
%
%   build_bleb_population_grid(mesh_dir, specs, 'CoreRangeNm', [10 120], ...
%       'BlebRangeNm', [10 120], 'NSizes', 22, 'OverlapFraction', 0.6)
%
%   mesh_dir : output folder (created if missing)
%   specs    : make_imaging_specs() output (uses specs.gap)
%
% Name-value options:
%   'CoreRangeNm'     : [min max] core diameter range (nm), default [10, 120].
%   'BlebRangeNm'     : [min max] bleb diameter range (nm), default [10, 120].
%   'NSizes'          : number of sizes per dimension, evenly spaced (linspace)
%                       across each range, default 22. Total particles = NSizes^2
%                       (default 22^2 = 484).
%   'OverlapFraction' : neck_overlap_frac, FIXED (not randomized) for every
%                       particle -- neck_overlap = OverlapFraction * min(d_core,
%                       d_bleb) / 2, exactly the live script's convention (there it
%                       was drawn uniformly per particle; here it's one fixed value
%                       so d_core/d_bleb are the only things varying). Default 0.6.
%   'ThetaTiltDeg'    : fixed tilt for every particle, default 90 -- "flat on the
%                       glass": the bleb-core neck axis lies entirely in the imaging
%                       (x-y) plane, so the bleb and core sit at the same z, rather
%                       than one stacked above the other. (0 = bleb directly above
%                       the core, along +z.)
%   'PhiAzimuthDeg'   : fixed azimuth for every particle, default 0. Can be changed
%                       per-particle AFTER the fact with adjust_bleb_azimuth.m --
%                       no need to rebuild the mesh for that (rebuilding IS needed
%                       to change ThetaTiltDeg; see that function's header for why).
%   'Step'            : voxelization step (nm), default 1 -- see build_fused_bleb_mesh.
%   'RadBound'        : target surface-triangle edge length (nm), default 6 -- but
%                       see 'MinRadBoundNm'/'RadBoundDivisor' below: this is a CEILING,
%                       used as-is only where the smaller of d_core/d_bleb is big
%                       enough for it, and automatically tightened for small sizes.
%   'MinRadBoundNm'   : floor on the adaptive radbound, default 3 (was 1.5 in an
%                       earlier version -- see the warning below for why that was
%                       too aggressive). Every grid cell actually uses radbound_eff =
%                       min(RadBound, max(MinRadBoundNm, min(d_core,d_bleb)/
%                       RadBoundDivisor)) -- a fixed 6 nm target cannot resolve the
%                       thin aqueous cap on a ~10-20 nm bleb (or core) as a separate
%                       labeled region from the lipid, and build_fused_bleb_mesh
%                       fails with "Only one labeled region was meshed" even though
%                       the two spheres genuinely do overlap (NOT true engulfment --
%                       just too coarse a target triangle size for that small a
%                       feature). Tightening radbound ONLY for the small end keeps
%                       runtime unchanged for most of the grid (cells with
%                       min(d_core,d_bleb) >= RadBound*RadBoundDivisor still get the
%                       plain 6 nm default) while making the small corner of the
%                       grid resolvable -- at the cost of more surface triangles
%                       (roughly (RadBound/radbound_eff)^2 more) and a correspondingly
%                       slower mesh + BEM solve for just those cells.
%                       WARNING: this only controls the MESH build here (this
%                       function never solves anything) -- but vol2mesh/cgalmesh
%                       takes one single global radbound for the WHOLE particle, so
%                       tightening it to resolve a small bleb ALSO refines the much
%                       larger core's surface, inflating TOTAL element count a lot
%                       more than the bleb alone would suggest (e.g. radbound 6->1.5
%                       nm on a ~100 nm core alone is a ~16x increase in core surface
%                       triangles). A cell meshed this fine can have tens of
%                       thousands of boundary elements, and nanobem's BEM solve
%                       builds an N-by-N array before it can even start solving --
%                       at N ~ 33000 that's over 8 GB and MATLAB errors out of
%                       BoundaryElement/touching (a repmat/ndgrid inside it) rather
%                       than solving. This function only BUILDS meshes, so it won't
%                       hit that error itself, but whatever solves these meshes
%                       later will. render_bleb_size_gradient.m guards against this
%                       with a post-build 'MaxElements' check (skips, doesn't crash)
%                       -- if you solve these meshes any other way, add the same
%                       guard, or just don't push MinRadBoundNm much below the
%                       default 3 without checking the resulting element counts.
%                       If a cell still fails at the floor, lower MinRadBoundNm
%                       further (slower, and more likely to hit the element-count
%                       wall above) or exclude that size from the sweep.
%   'RadBoundDivisor' : default 5 -- radbound_eff targets roughly this many triangle
%                       edges across the smaller of d_core/d_bleb.
%   'Prefix'          : filename prefix, default 'blebgrid'. Files are saved as
%                       <prefix>_<core_index>_<bleb_index>_mesh.mat with 2-digit
%                       zero-padded grid indices (e.g. 'blebgrid_01_01_mesh.mat'
%                       .. 'blebgrid_22_22_mesh.mat'), so they sort/glob cleanly and
%                       each filename alone says which grid cell it is.
%
% Each saved mesh.params also carries overlap_fraction, core_index, bleb_index,
% core_size_nm, bleb_size_nm, n_sizes, so any file alone documents its place in the
% grid without needing to re-derive it from the filename.
%
% Resumable: an existing file is skipped without recomputation (build_fused_bleb_mesh
% -- the voxelize + vol2mesh step -- is what's slow here; there is nothing random to
% re-draw on a rerun, unlike the seeded-RNG population builders).
%
% A geometrically-invalid cell (bleb too thin to resolve at this Step/OverlapFraction
% -- see build_fused_bleb_mesh's error) is logged and skipped, not retried: with
% every parameter fixed per cell there is nothing stochastic left to redraw, so a
% retry would just fail the same way. If a run still reports "Only one labeled
% region was meshed" failures after the adaptive radbound above, lower
% 'MinRadBoundNm' (slower) or 'OverlapFraction', or shrink 'Step'; rerun (already
% -saved cells are skipped).
%
% Requires iso2mesh on the path (same as build_fused_bleb_mesh.m).

    p = inputParser;
    addParameter(p, 'CoreRangeNm', [10, 120]);
    addParameter(p, 'BlebRangeNm', [10, 120]);
    addParameter(p, 'NSizes', 22);
    addParameter(p, 'OverlapFraction', 0.6);
    addParameter(p, 'ThetaTiltDeg', 90);
    addParameter(p, 'PhiAzimuthDeg', 0);
    addParameter(p, 'Step', 1);
    addParameter(p, 'RadBound', 6);
    addParameter(p, 'MinRadBoundNm', 3);
    addParameter(p, 'RadBoundDivisor', 5);
    addParameter(p, 'Prefix', 'blebgrid');
    parse(p, varargin{:});
    opt = p.Results;

    if ~exist(mesh_dir, 'dir'), mkdir(mesh_dir); end

    n = opt.NSizes;
    core_sizes = linspace(opt.CoreRangeNm(1), opt.CoreRangeNm(2), n);
    bleb_sizes = linspace(opt.BlebRangeNm(1), opt.BlebRangeNm(2), n);

    fprintf('build_bleb_population_grid: %d core sizes x %d bleb sizes = %d particles\n', ...
        n, n, n^2);
    fprintf('  core sizes (nm): %s\n', mat2str(round(core_sizes, 2)));
    fprintf('  bleb sizes (nm): %s\n', mat2str(round(bleb_sizes, 2)));
    fprintf('  OverlapFraction = %.3g, ThetaTiltDeg = %.3g, PhiAzimuthDeg = %.3g\n', ...
        opt.OverlapFraction, opt.ThetaTiltDeg, opt.PhiAzimuthDeg);

    n_ok = 0; n_skip = 0; n_fail = 0;
    failed = {};

    for ic = 1:n
        for ib = 1:n
            d_core = core_sizes(ic);
            d_bleb = bleb_sizes(ib);
            fn = fullfile(mesh_dir, sprintf('%s_%02d_%02d_mesh.mat', opt.Prefix, ic, ib));

            if exist(fn, 'file')
                n_skip = n_skip + 1;
                continue;
            end

            neck_overlap = opt.OverlapFraction * min(d_core, d_bleb) / 2;
            if d_bleb <= neck_overlap + 4 * opt.Step
                warning('build_bleb_population_grid:tooThin', ...
                    ['[%d,%d] d_core=%.2f d_bleb=%.2f: bleb too thin at this overlap/step ' ...
                     '(neck_overlap=%.2f) -- skipped.'], ic, ib, d_core, d_bleb, neck_overlap);
                n_fail = n_fail + 1;
                failed{end+1} = sprintf('[%d,%d] d_core=%.2f d_bleb=%.2f: too thin (neck_overlap=%.2f)', ...
                    ic, ib, d_core, d_bleb, neck_overlap); %#ok<AGROW>
                continue;
            end

            radbound_eff = min(opt.RadBound, ...
                max(opt.MinRadBoundNm, min(d_core, d_bleb) / opt.RadBoundDivisor));

            try
                mesh = build_fused_bleb_mesh(d_core, d_bleb, neck_overlap, opt.Step, ...
                    radbound_eff, specs.gap, opt.ThetaTiltDeg, opt.PhiAzimuthDeg);
            catch ME
                warning('build_bleb_population_grid:buildFailed', ...
                    '[%d,%d] d_core=%.2f d_bleb=%.2f failed: %s', ic, ib, d_core, d_bleb, ME.message);
                n_fail = n_fail + 1;
                failed{end+1} = sprintf('[%d,%d] d_core=%.2f d_bleb=%.2f: %s', ...
                    ic, ib, d_core, d_bleb, ME.message); %#ok<AGROW>
                continue;
            end

            mesh.params.overlap_fraction = opt.OverlapFraction;
            mesh.params.core_index       = ic;
            mesh.params.bleb_index       = ib;
            mesh.params.core_size_nm     = d_core;
            mesh.params.bleb_size_nm     = d_bleb;
            mesh.params.n_sizes          = n;

            save(fn, 'mesh');
            n_ok = n_ok + 1;
            fprintf('saved %s (d_core=%.2f, d_bleb=%.2f)\n', fn, d_core, d_bleb);
        end
    end

    fprintf('\nbuild_bleb_population_grid done: %d saved, %d already existed, %d failed/skipped.\n', ...
        n_ok, n_skip, n_fail);
    if n_fail > 0
        fprintf('failed/skipped grid cells:\n');
        for i = 1:numel(failed)
            fprintf('  %s\n', failed{i});
        end
    end
end
