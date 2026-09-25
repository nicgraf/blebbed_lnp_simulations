function out = render_mesh_size_gallery(grid_dir, varargin)
% RENDER_MESH_SIZE_GALLERY  Sample particles across a build_bleb_population_grid.m
% library and draw their actual 3D mesh surfaces side by side, all in ONE SHARED
% physical scale, to show how much particle size varies across the library at a
% glance (e.g. for a slide). No BEM solving involved -- this only loads
% mesh.node/f_lipid_out/f_aqueous_out/f_interface and draws them; it's fast even
% for a library that hasn't been solved at all yet.
%
%   render_mesh_size_gallery(grid_dir)
%       -> default 22x22 grid (484 particles total, blebgrid_<core>_<bleb>_mesh.mat),
%          targets every 15th particle by flat index (15, 30, 45, ..., 480 -- 32
%          targets, floor(484/15) = 32), one 3D-rendered panel each, saved as a PNG.
%
%   render_mesh_size_gallery(grid_dir, 'Step', 20)   % every 20th -> 24 targets
%
% Indexing: particles are numbered 1..NSizes^2 in the same nested-loop order
% build_bleb_population_grid.m saves them in (core index outer, bleb index inner):
% flat index i -> core_index = ceil(i/NSizes), bleb_index = i - (core_index-1)*NSizes.
% So particle 1 is <Prefix>_01_01, particle 22 is <Prefix>_01_22, particle 23 is
% <Prefix>_02_01, and so on.
%
% If a target index's mesh file doesn't exist yet (the grid isn't fully built),
% this does NOT leave a gap or build it -- it steps forward one particle at a time
% (index+1, index+2, ...) until it finds one that DOES exist, and plots that one
% instead. The search only ever moves forward and never reuses a particle already
% shown by an earlier target, so every panel that's supposed to render, does --
% you always get up to nSample real particles (fewer only if the library runs out
% of remaining particles near the very end). out.idx reports which particle
% actually ended up in each panel, which can differ from the original target grid
% (Step, 2*Step, ...) if any nudging happened.
%
% grid_dir : the folder build_bleb_population_grid.m saved the library into.
%
% Name-value options:
%   'Prefix'      : filename prefix, default 'blebgrid' -- must match the grid's Prefix.
%   'NSizes'      : sizes per dimension, default 22 -- must match the grid's NSizes
%                   (total particles sampled from = NSizes^2).
%   'CoreRangeNm' , 'BlebRangeNm' : [min max] ranges (nm), default [10, 120] each --
%                   must match the grid's, only used to print/label actual sizes.
%   'Step'        : target every this-many-th particle by flat index, default 15.
%                   Targets are Step, 2*Step, ..., floor(NSizes^2/Step)*Step -- e.g.
%                   the default gives 32 targets (floor(484/15) = 32) at 15, 30,
%                   ..., 480. A missing target is nudged forward to the next
%                   existing particle (see above) -- so this is the requested
%                   sampling density, not a guaranteed exact index.
%   'FaceColorLipid'   : RGB, default [0.95 0.85 0.65] (warm) for the lipid core surface.
%   'FaceColorAqueous' : RGB, default [0.65 0.80 0.95] (cool) for the aqueous bleb
%                        surface AND the lipid/aqueous interface (same material side
%                        visually -- it's the boundary of the same aqueous pocket).
%   'GridSize'    : [nrows ncols], default [] = auto, wider than tall (~1.3x aspect),
%                   sized to just fit the sampled count.
%   'ViewAngle'   : [az el] for every panel's 3D view, default [-35 20].
%   'SaveFigure'  : default true.
%   'FigureFile'  : default '' = fullfile(grid_dir, 'mesh_size_gallery.png'), written
%                   via exportgraphics at 300 dpi.
%
% All panels share ONE physical scale (the same x/y/z axis limits, sized to the
% largest sampled particle's bounding box) -- this is what makes the size
% differences actually visible: a panel is not auto-cropped to its own particle,
% so a 20 nm particle looks small and a 200 nm one fills the frame, exactly as they
% would side by side in reality. (Every particle's lipid core is centered at
% x=y=0 by construction -- see build_fused_bleb_mesh.m / adjust_bleb_azimuth.m's
% header -- so a single shared frame centered at the origin works for all of them
% without per-panel re-centering.)
%
% A target is only left unfilled (panel shows "missing") if the forward search
% runs off the end of the library (index > NSizes^2) without finding another
% existing file -- meaning the grid is missing a whole trailing stretch, not just
% one cell. Run build_bleb_population_grid.m first if the grid isn't fully
% populated yet.
%
% out fields:
%   idx                        : the ACTUAL flat index shown in each panel (1-based,
%                                into 1..NSizes^2) -- may differ from the original
%                                Step-spaced targets where a target was nudged forward
%   coreIndex, blebIndex        : grid (ic, ib) for each panel's actual index
%   coreSizeNm, blebSizeNm      : actual sizes (nm) for each panel's actual index
%   ok                          : logical, false only where no existing particle
%                                was found for that panel at all (ran off the end)
%   fig                         : the figure handle
%   figureFile                  : path written, or '' if SaveFigure is false

    p = inputParser;
    addParameter(p, 'Prefix', 'blebgrid');
    addParameter(p, 'NSizes', 22);
    addParameter(p, 'CoreRangeNm', [10, 120]);
    addParameter(p, 'BlebRangeNm', [10, 120]);
    addParameter(p, 'Step', 15);
    addParameter(p, 'FaceColorLipid', [0.95 0.85 0.65]);
    addParameter(p, 'FaceColorAqueous', [0.65 0.80 0.95]);
    addParameter(p, 'GridSize', []);
    addParameter(p, 'ViewAngle', [-35 20]);
    addParameter(p, 'SaveFigure', true);
    addParameter(p, 'FigureFile', '');
    parse(p, varargin{:});
    opt = p.Results;

    assert(exist(grid_dir, 'dir') == 7, 'render_mesh_size_gallery:noGridDir', ...
        '%s does not exist.', grid_dir);

    nTotal = opt.NSizes^2;
    nSample = floor(nTotal / opt.Step);
    assert(nSample >= 1, 'render_mesh_size_gallery:stepTooLarge', ...
        'Step=%d leaves no samples out of %d particles.', opt.Step, nTotal);
    targets = (1:nSample) * opt.Step;

    core_sizes = linspace(opt.CoreRangeNm(1), opt.CoreRangeNm(2), opt.NSizes);
    bleb_sizes = linspace(opt.BlebRangeNm(1), opt.BlebRangeNm(2), opt.NSizes);

    fprintf(['render_mesh_size_gallery: %d particles total, every %d -> %d targets ' ...
        '(%s) -- a missing target is nudged forward to the next existing particle\n'], ...
        nTotal, opt.Step, nSample, mat2str(targets));

    % ---- resolve each target to an actual (existing) particle, monotonically
    % forward, no repeats -- then load it and track the shared bounding box ----
    idx        = nan(1, nSample);
    coreIndex  = nan(1, nSample);
    blebIndex  = nan(1, nSample);
    coreSizeNm = nan(1, nSample);
    blebSizeNm = nan(1, nSample);
    meshList   = cell(1, nSample);
    ok         = false(1, nSample);
    globalMin  = [Inf Inf Inf];
    globalMax  = [-Inf -Inf -Inf];

    cursor = 0;   % last particle actually used -- never reuse or go backward
    for k = 1:nSample
        cand = max(targets(k), cursor + 1);
        found = false;
        while cand <= nTotal
            ic = ceil(cand / opt.NSizes);
            ib = cand - (ic - 1) * opt.NSizes;
            fn = fullfile(grid_dir, sprintf('%s_%02d_%02d_mesh.mat', opt.Prefix, ic, ib));
            if exist(fn, 'file')
                found = true;
                break;
            end
            cand = cand + 1;
        end

        if ~found
            warning('render_mesh_size_gallery:ranOut', ...
                ['no existing particle found at or after target %d for panel %d/%d -- ' ...
                 'the library is missing a trailing stretch; leaving this panel blank.'], ...
                targets(k), k, nSample);
            continue;
        end

        if cand > targets(k)
            fprintf('  target %d missing -- nudged forward to particle %d\n', targets(k), cand);
        end

        M = load(fn, 'mesh');
        meshList{k} = M.mesh;
        ok(k) = true;
        idx(k) = cand;
        coreIndex(k) = ic;
        blebIndex(k) = ib;
        coreSizeNm(k) = core_sizes(ic);
        blebSizeNm(k) = bleb_sizes(ib);
        cursor = cand;
        globalMin = min(globalMin, min(M.mesh.node, [], 1));
        globalMax = max(globalMax, max(M.mesh.node, [], 1));
        fprintf('  [%d/%d] particle %d: core %.1f nm / bleb %.1f nm (%s)\n', ...
            k, nSample, cand, coreSizeNm(k), blebSizeNm(k), fn);
    end

    if ~any(ok)
        warning('render_mesh_size_gallery:noneFound', ...
            'no existing mesh files were found for any of the %d targets -- nothing to plot.', ...
            nSample);
    end

    % ---- one shared physical frame for every panel (centered at x=y=0) ----
    pad = 1.08;
    M_xy = pad * max([abs(globalMin(1:2)), abs(globalMax(1:2))]);
    if ~isfinite(M_xy), M_xy = 1; end
    zTop = pad * globalMax(3);
    if ~isfinite(zTop), zTop = 1; end

    % ---- layout ----
    if isempty(opt.GridSize)
        ncols = max(1, round(sqrt(nSample * 1.3)));
        ncols = min(ncols, nSample);
        nrows = ceil(nSample / ncols);
    else
        nrows = opt.GridSize(1);
        ncols = opt.GridSize(2);
        assert(nrows * ncols >= nSample, 'render_mesh_size_gallery:gridTooSmall', ...
            '''GridSize'' [%d %d] has only %d tiles for %d panels.', ...
            nrows, ncols, nrows * ncols, nSample);
    end

    fig = figure('Color', 'w', 'Position', [50 50 200 * ncols 200 * nrows]);
    tl = tiledlayout(fig, nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');

    for k = 1:nSample
        ax = nexttile(tl);
        if ~ok(k)
            axis(ax, 'off');
            text(ax, 0.5, 0.5, 'missing', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', [0.5 0.5 0.5]);
            continue;
        end

        mesh = meshList{k};
        hold(ax, 'on');
        if ~isempty(mesh.f_lipid_out)
            trisurf(double(mesh.f_lipid_out), mesh.node(:,1), mesh.node(:,2), mesh.node(:,3), ...
                'Parent', ax, 'FaceColor', opt.FaceColorLipid, 'EdgeColor', 'none', ...
                'FaceLighting', 'gouraud');
        end
        if ~isempty(mesh.f_aqueous_out)
            trisurf(double(mesh.f_aqueous_out), mesh.node(:,1), mesh.node(:,2), mesh.node(:,3), ...
                'Parent', ax, 'FaceColor', opt.FaceColorAqueous, 'EdgeColor', 'none', ...
                'FaceLighting', 'gouraud');
        end
        if ~isempty(mesh.f_interface)
            trisurf(double(mesh.f_interface), mesh.node(:,1), mesh.node(:,2), mesh.node(:,3), ...
                'Parent', ax, 'FaceColor', opt.FaceColorAqueous, 'EdgeColor', 'none', ...
                'FaceLighting', 'gouraud');
        end
        hold(ax, 'off');

        axis(ax, 'equal');
        xlim(ax, [-M_xy, M_xy]);
        ylim(ax, [-M_xy, M_xy]);
        zlim(ax, [0, zTop]);
        axis(ax, 'off');
        view(ax, opt.ViewAngle);
        camlight(ax, 'headlight');
        lighting(ax, 'gouraud');
        title(ax, sprintf('%.0f / %.0f nm', coreSizeNm(k), blebSizeNm(k)), ...
            'FontSize', 8, 'FontWeight', 'normal');
    end
    title(tl, sprintf('Mesh library sample (every %d of %d particles -- core / bleb, nm)', ...
        opt.Step, nTotal));

    figureFile = '';
    if opt.SaveFigure
        figureFile = opt.FigureFile;
        if isempty(figureFile)
            figureFile = fullfile(grid_dir, 'mesh_size_gallery.png');
        end
        exportgraphics(fig, figureFile, 'Resolution', 300);
        fprintf('saved %s\n', figureFile);
    end

    out = struct('idx', idx, 'coreIndex', coreIndex, 'blebIndex', blebIndex, ...
        'coreSizeNm', coreSizeNm, 'blebSizeNm', blebSizeNm, 'ok', ok, ...
        'fig', fig, 'figureFile', figureFile);
end
