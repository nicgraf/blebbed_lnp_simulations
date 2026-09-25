function out = render_bleb_size_gradient(grid_dir, specs, varargin)
% RENDER_BLEB_SIZE_GRADIENT  One row of a build_bleb_population_grid.m population --
% a fixed core size, every bleb size on the grid -- laid out as in-plane iSCAT
% images in one grid figure, e.g. to show "how does the bleb signature change with
% bleb size" at a glance, for a slide.
%
% This function does NOT build meshes. It expects the grid cells for the requested
% core row (e.g. blebgrid_18_01_mesh.mat .. blebgrid_18_22_mesh.mat) to already
% exist in grid_dir, from build_bleb_population_grid.m. A cell whose mesh file is
% missing is skipped (with a clear reason), not built -- run
% build_bleb_population_grid.m first if you need to fill it in.
%
% Two strictly separate passes, each keeping only one particle's data in memory
% at a time:
%   1. SOLVE PASS: for bleb index 1..NSizes, load that cell's mesh (not the whole
%      row), rebuild its BEM boundary, solve, and save the solution back into that
%      cell's .mat file IMMEDIATELY -- then everything about that particle (mesh,
%      boundary elements, solution) is cleared before moving to the next index. A
%      cell that's already solved (bem_sol present in the file) is skipped without
%      loading it. Nothing from this pass is kept in memory once it's done.
%   2. IMAGING PASS: for bleb index 1..NSizes, reload just that cell's saved
%      solution, propagate the in-plane image at ZFocusNm, store the resulting
%      (small) contrast image, then discard the solution before the next index.
% Only after both passes finish does the figure get built from the collected
% images. This means: (a) at no point does the function hold more than one
% particle's mesh/solution in memory, however many sizes are in the row, and
% (b) the two passes are independently rerunnable -- e.g. call again with a
% different ZFocusNm and it reuses every already-solved cell, redoing only the
% (cheap) imaging pass.
%
%   render_bleb_size_gradient(grid_dir, specs)
%       -> grid_dir must be the SAME folder you passed to build_bleb_population_grid.
%          Default CoreSizeNm = 100 nm resolves to the NEAREST core size actually on
%          the grid (with the default 22-point 10-120 nm grid, that's core index 18,
%          ~99.05 nm -- printed, not silently substituted) and sweeps that row's 22
%          bleb sizes, z focus = 50 nm -- 22 panels in a 4x6 grid, one shared color
%          scale, saved as a PNG next to the meshes.
%
%   render_bleb_size_gradient(grid_dir, specs, 'CoreIndex', 18)   % same row, by index
%   render_bleb_size_gradient(grid_dir, specs, 'ZFocusNm', 0, 'FigureFile', 'z0.png')
%
%   grid_dir : the SAME folder used with build_bleb_population_grid.m (e.g. your
%              'testing_grid' folder). Must already exist and already contain the
%              requested row's mesh files.
%   specs    : make_imaging_specs() output.
%
% Name-value options (defaults match build_bleb_population_grid.m's -- keep
% 'Prefix'/'CoreRangeNm'/'BlebRangeNm'/'NSizes' the same as whatever you actually
% built the grid with, or the row/column indices won't line up with the right sizes):
%   'CoreSizeNm'      : the core size (nm) you want, default 100 -- resolved to the
%                       NEAREST index on the grid (core_sizes = linspace(CoreRangeNm(1),
%                       CoreRangeNm(2), NSizes)); ignored if 'CoreIndex' is given.
%   'CoreIndex'       : grid row index directly (1..NSizes), default [] = resolve
%                       from CoreSizeNm instead. Use this if you already know the
%                       index (e.g. from a previous run's printed resolution).
%   'CoreRangeNm'     : [min max] core diameter range (nm), default [10, 120] --
%                       must match the grid's CoreRangeNm (only used to resolve
%                       CoreSizeNm/print sizes -- filenames are by index).
%   'BlebRangeNm'     : [min max] bleb diameter range (nm), default [10, 120] --
%                       must match the grid's BlebRangeNm.
%   'NSizes'          : sizes per dimension, default 22 -- must match the grid's NSizes.
%   'Prefix'          : filename prefix, default 'blebgrid' -- must match the grid's
%                       Prefix. Files are read as <Prefix>_<CoreIndex>_<ib>_mesh.mat,
%                       exactly build_bleb_population_grid.m's convention.
%   'MaxElements'     : safety cap on the number of BEM boundary elements a cell is
%                       allowed to need, default 12000. This exists because
%                       nanobem's BEM solve builds an N-by-N working array (in
%                       BoundaryElement/touching, before the actual solve), and that
%                       array's memory is EXACTLY N^2 * ~7.5 bytes: at N=12000 that's
%                       ~1.1 GB (safe headroom below MATLAB's default 8 GB array-size
%                       preference, which itself exists to stop MATLAB going
%                       unresponsive). A cell meshed too fine (e.g. from an old,
%                       more aggressive MinRadBoundNm in build_bleb_population_grid.m)
%                       can have tens of thousands of elements -- N > 33000 needs a
%                       >8 GB array and MATLAB errors out of BoundaryElement/touching
%                       rather than solving. A cell whose ACTUAL element count
%                       (checked right after rebuild_bleb_tau, before the solve)
%                       exceeds MaxElements is SKIPPED in the solve pass (logged, not
%                       crashed) -- see out.skipReason and the placeholder tiles in
%                       the figure. Raise this only if you know your machine has the
%                       RAM; a crashed/unresponsive MATLAB is worse than a skipped
%                       panel.
%   'ZFocusNm'        : focus height (nm) at which each in-plane image is taken,
%                       default 50.
%   'GridSize'        : [nrows ncols] for the panel layout, default [] = auto
%                       (ceil(NSizes/6)-by-6, i.e. [4 6] for the default 22 sizes).
%   'SaveFigure'       : default true.
%   'FigureFile'       : path to save the figure to (.png/.pdf/...), default '' =
%                       auto: fullfile(grid_dir, 'bleb_size_gradient_core<CoreIndex>.png').
%                       Written via exportgraphics at 300 dpi -- crisp when pasted
%                       into a slide at typical sizes. Ignored if SaveFigure is false.
%
% A cell that fails for ANY reason (mesh file missing, too many elements for
% MaxElements, or anything else) is logged with a warning and SKIPPED -- it does
% not abort the run. Its panel in the figure is a light gray placeholder with a
% short wrapped reason instead of an image, so a partial result is still
% immediately readable rather than all-or-nothing. See out.skipReason (one cell
% per bleb size, '' if OK) and the printed summary at the end of the run.
%
% Every panel with real data shares ONE color scale (cube-root-compressed, like
% set_cbrt_colorbar.m, with a single colorbar for the whole grid) -- so panels are
% directly comparable to each other and the trend across bleb sizes isn't an
% artifact of each panel re-normalizing to its own contrast range.
%
% Runtime: solving is done inline per particle (order = specs.bem_order stratified
% solve each) -- the slow part, likely minutes per particle -- and does NOT use
% solve_bleb_population.m (which has an unrelated hardcoded `for f = 11:50` loop
% range that would silently under-solve a batch this size, AND no MaxElements-style
% guard against this same out-of-memory failure; worth fixing separately, flagged,
% not touched here). If the whole row is already solved, the solve pass is a fast
% no-op (just a variable-list check per file) and only the imaging pass runs.
%
% out fields:
%   d_bleb                    : the NSizes swept bleb diameters (nm)
%   contrast                   : NSizes-by-ny-by-nx contrast image stack (%); a
%                                skipped row's slice is left as zeros -- check
%                                skipReason, don't infer failure from an all-zero slice
%   skipReason                 : 1-by-NSizes cell array, '' for a panel that rendered,
%                                else the reason it was skipped
%   cmax                       : the shared |contrast| color limit used across all
%                                panels that rendered
%   core_index, core_size_nm, z_focus_nm : as resolved/used
%   fig, tl                    : the figure and tiledlayout handles
%   figureFile                 : path written, or '' if SaveFigure is false

    p = inputParser;
    addParameter(p, 'CoreSizeNm', 100);
    addParameter(p, 'CoreIndex', []);
    addParameter(p, 'CoreRangeNm', [10, 120]);
    addParameter(p, 'BlebRangeNm', [10, 120]);
    addParameter(p, 'NSizes', 22);
    addParameter(p, 'Prefix', 'blebgrid');
    addParameter(p, 'MaxElements', 12000);
    addParameter(p, 'ZFocusNm', 50);
    addParameter(p, 'GridSize', []);
    addParameter(p, 'SaveFigure', true);
    addParameter(p, 'FigureFile', '');
    parse(p, varargin{:});
    opt = p.Results;

    assert(exist(grid_dir, 'dir') == 7, 'render_bleb_size_gradient:noGridDir', ...
        ['%s does not exist -- this reuses an existing build_bleb_population_grid ' ...
         'folder, it does not create a new one. Build the grid first, or point ' ...
         '''grid_dir'' at the folder you already built.'], grid_dir);
    if ~isfield(specs, 'tau_method') || isempty(specs.tau_method), specs.tau_method = 'native'; end
    sys = build_imaging_system(specs);

    n = opt.NSizes;
    core_sizes = linspace(opt.CoreRangeNm(1), opt.CoreRangeNm(2), n);
    bleb_sizes = linspace(opt.BlebRangeNm(1), opt.BlebRangeNm(2), n);

    if isempty(opt.CoreIndex)
        [~, coreIndex] = min(abs(core_sizes - opt.CoreSizeNm));
        d_core = core_sizes(coreIndex);
        if abs(d_core - opt.CoreSizeNm) > 1e-9
            fprintf(['render_bleb_size_gradient: CoreSizeNm=%.2f nm is not exactly on the ' ...
                'grid -- using the nearest grid point, core index %d (%.3f nm).\n'], ...
                opt.CoreSizeNm, coreIndex, d_core);
        end
    else
        coreIndex = opt.CoreIndex;
        d_core = core_sizes(coreIndex);
    end

    fns = cell(1, n);
    for ib = 1:n
        fns{ib} = fullfile(grid_dir, sprintf('%s_%02d_%02d_mesh.mat', opt.Prefix, coreIndex, ib));
    end
    skipReason = cell(1, n);

    fprintf(['render_bleb_size_gradient: core index %d/%d (%.2f nm), %d bleb sizes ' ...
        '%.1f-%.1f nm, z = %.0f nm, %s/%s_%02d_*_mesh.mat\n'], ...
        coreIndex, n, d_core, n, bleb_sizes(1), bleb_sizes(end), opt.ZFocusNm, ...
        grid_dir, opt.Prefix, coreIndex);

    % ============================================================
    % PASS 1 -- solve every missing cell, save immediately, discard.
    % Holds at most ONE particle's mesh/tau/solution in memory at a time.
    % ============================================================
    fprintf('-- solve pass --\n');
    for ib = 1:n
        fn = fns{ib};
        d_bleb = bleb_sizes(ib);

        if ~exist(fn, 'file')
            skipReason{ib} = sprintf(['mesh file not found (%s) -- this function does not ' ...
                'build meshes; run build_bleb_population_grid.m first.'], fn);
            warning('render_bleb_size_gradient:meshMissing', '[%d/%d] %s', ib, n, skipReason{ib});
            continue;
        end

        if ismember('bem_sol', who('-file', fn))
            fprintf('  [%d/%d] d_bleb=%.2f nm already solved, skipping\n', ib, n, d_bleb);
            continue;
        end

        try
            M = load(fn, 'mesh');
            tau = rebuild_bleb_tau(M.mesh, sys.mat_set, sys.idx_water, sys.idx_lipid, ...
                sys.idx_aqueous, sys.tau_method);
            clear M

            if numel(tau) > opt.MaxElements
                error('render_bleb_size_gradient:tooManyElements', ...
                    ['%d boundary elements > MaxElements=%d -- the BEM solve''s N-by-N ' ...
                     'working array would need roughly %.1f GB. Raise MaxElements only ' ...
                     'if you know you have the RAM; otherwise this size needs a coarser ' ...
                     'mesh (rebuild that cell with a larger RadBound/MinRadBoundNm) or ' ...
                     'should be dropped from the sweep.'], ...
                    numel(tau), opt.MaxElements, (double(numel(tau))^2 * 7.5) / 1e9);
            end

            qinc = sys.einc(tau, 'layer', sys.layer);
            bem  = stratified.bemsolver(tau, sys.layer, 'order', specs.bem_order);
            sol  = bem \ qinc;
            bem_sol = struct('e', sol.e, 'h', sol.h, 'k0', sol.k0, ...
                'nelem', numel(tau), 'pos_checksum', sum(vertcat(tau.pos), 'all'), ...
                'lambda', specs.lambda, 'pol', specs.pol, 'dir', specs.dir);
            save(fn, 'bem_sol', '-append');
            fprintf('  [%d/%d] d_bleb=%.2f nm solved + saved (%d elements)\n', ...
                ib, n, d_bleb, numel(tau));
            clear tau qinc bem sol bem_sol

        catch ME
            warning('render_bleb_size_gradient:cellSkipped', ...
                '[%d/%d] d_bleb=%.2f nm skipped: %s', ib, n, d_bleb, ME.message);
            skipReason{ib} = ME.message;
            clear M tau qinc bem sol
        end
    end

    % ============================================================
    % PASS 2 -- reload one saved solution at a time, image it, discard.
    % ============================================================
    fprintf('-- imaging pass (z = %.0f nm) --\n', opt.ZFocusNm);
    x = sys.x;
    nxpix = numel(x);
    contrastStack = zeros(n, nxpix, nxpix);

    for ib = 1:n
        if ~isempty(skipReason{ib}), continue; end
        fn = fns{ib};
        d_bleb = bleb_sizes(ib);
        try
            [~, sol] = load_bleb_solution(fn, sys);
            far = farfields(sol, sys.lens.dir);
            contrast = propagate_iscat_image(sol, opt.ZFocusNm, sys, far);
            contrastStack(ib, :, :) = contrast;
            fprintf('  [%d/%d] d_bleb=%.2f nm imaged\n', ib, n, d_bleb);
            clear sol far contrast
        catch ME
            warning('render_bleb_size_gradient:imagingSkipped', ...
                '[%d/%d] d_bleb=%.2f nm skipped during imaging: %s', ib, n, d_bleb, ME.message);
            skipReason{ib} = ME.message;
        end
    end

    ok = cellfun(@isempty, skipReason);
    nSkipped = sum(~ok);

    % ---- one shared color scale across every rendered panel ----
    if any(ok)
        cmax = max(abs(contrastStack(ok, :, :)), [], 'all');
    else
        cmax = 0;
    end
    if cmax == 0, cmax = 1; end

    % ---- grid figure ----
    if isempty(opt.GridSize)
        ncols = min(n, 6);
        nrows = ceil(n / ncols);
    else
        nrows = opt.GridSize(1);
        ncols = opt.GridSize(2);
        if nrows * ncols < n
            error('render_bleb_size_gradient:gridTooSmall', ...
                '''GridSize'' [%d %d] has only %d tiles for %d panels.', ...
                nrows, ncols, nrows * ncols, n);
        end
    end

    fig = figure('Color', 'w', 'Position', [100 100 220 * ncols 240 * nrows]);
    colormap(fig, 'gray');   % symmetric gray LUT: 0 contrast = mid-gray, shared by every tile
    tl = tiledlayout(fig, nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');
    lastAx = [];
    for ib = 1:n
        ax = nexttile(tl);
        if ok(ib)
            img = squeeze(contrastStack(ib, :, :));
            imagesc(ax, x, x, nthroot(img, 3), nthroot([-cmax, cmax], 3));
            axis(ax, 'image');
            lastAx = ax;
        else
            axis(ax, 'image');
            xlim(ax, [0 1]); ylim(ax, [0 1]);
            patch(ax, [0 1 1 0], [0 1 0 1], [0.92 0.92 0.9], 'EdgeColor', 'none');
            reasonShort = skipReason{ib};
            if numel(reasonShort) > 40, reasonShort = [reasonShort(1:37) '...']; end
            text(ax, 0.5, 0.5, sprintf('skipped\n%s', reasonShort), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', 7, 'Color', [0.4 0.4 0.38], 'Interpreter', 'none');
        end
        ax.XTick = [];
        ax.YTick = [];
        box(ax, 'on');
        ax.XColor = [0.6 0.6 0.6];
        ax.YColor = [0.6 0.6 0.6];
        title(ax, sprintf('%.0f nm', bleb_sizes(ib)), 'FontWeight', 'normal', 'FontSize', 9);
    end
    title(tl, sprintf('Bleb size gradient (core = %.2f nm, z = %.0f nm)', d_core, opt.ZFocusNm));

    % ---- one shared colorbar for the whole grid, cube-root ticks (like set_cbrt_colorbar.m) ----
    if ~isempty(lastAx)
        cb = colorbar(lastAx);
        cb.Layout.Tile = 'east';
        decade_ticks = [1 2 5];
        all_ticks = [];
        for e = 0:4, all_ticks = [all_ticks, decade_ticks * 10^e]; end %#ok<AGROW>
        all_ticks = unique(all_ticks);
        all_ticks = all_ticks(all_ticks <= cmax);
        if isempty(all_ticks), all_ticks = cmax; end
        tick_vals = [-fliplr(all_ticks), 0, all_ticks];
        cb.Ticks = nthroot(tick_vals, 3);
        cb.TickLabels = arrayfun(@(v) sprintf('%g%%', v), tick_vals, 'UniformOutput', false);
        cb.Label.String = 'iSCAT contrast (%)';
    end

    figureFile = '';
    if opt.SaveFigure
        figureFile = opt.FigureFile;
        if isempty(figureFile)
            figureFile = fullfile(grid_dir, sprintf('bleb_size_gradient_core%02d.png', coreIndex));
        end
        exportgraphics(fig, figureFile, 'Resolution', 300);
        fprintf('saved %s\n', figureFile);
    end

    fprintf('render_bleb_size_gradient done: %d/%d panels rendered', n - nSkipped, n);
    if nSkipped > 0
        fprintf(', %d skipped:\n', nSkipped);
        for ib = 1:n
            if ~ok(ib)
                fprintf('  index %d (%.2f nm): %s\n', ib, bleb_sizes(ib), skipReason{ib});
            end
        end
    else
        fprintf('.\n');
    end

    out = struct('d_bleb', bleb_sizes, 'contrast', contrastStack, 'skipReason', {skipReason}, ...
        'cmax', cmax, 'core_index', coreIndex, 'core_size_nm', d_core, 'z_focus_nm', opt.ZFocusNm, ...
        'fig', fig, 'tl', tl, 'figureFile', figureFile);
end
