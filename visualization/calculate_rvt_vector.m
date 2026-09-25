function vec = calculate_rvt_vector(rvt_out, varargin)
% CALCULATE_RVT_VECTOR  Trace the RVT-map peak (the "most symmetric point") across
% EVERY z plane of an already-computed RVT z-stack: at each z, lightly Gaussian-smooth
% the RVT map, find the (x,y) position and value of its maximum, and return/plot that
% trace as-is -- with no assumption that it moves in a straight line. (An earlier
% version of this function collapsed the trace to a single start->end "net vector" and
% plotted only that; this version plots the real point-by-point path instead, since the
% actual motion can double back, jump between competing symmetric features, sit near the
% FOV edge, etc. -- collapsing it to one arrow hid exactly the behavior worth looking at.)
%
%   vec = calculate_rvt_vector(rvt_out)
%   vec = calculate_rvt_vector(rvt_out, 'Display', true)
%
%   % overlay a sphere's and a bleb's traces on the same 4 panels for direct comparison:
%   vecS = calculate_rvt_vector(outS, 'Display', true);
%   vecB = calculate_rvt_vector(outB, 'Display', true, 'Axes', vecS.ax, ...
%              'Color', '#eb6834', 'DisplayName', 'bleb 001');
%
%   rvt_out : the struct returned by propagate_iscat_rvt_zstack, called with its default
%             'StoreMaps', true (so rvt_out.rvtMaps and rvt_out.x are present). The peak
%             position is re-derived from the raw maps on every call here, independent
%             of whatever 'CenterMode' the earlier propagate_iscat_rvt_zstack call used.
%             rvt_out.rvtMaps itself is never modified -- the smoothing below (see
%             'SmoothSigmaPx') is applied to a local copy of each map, used only to find
%             that map's peak, so a second call (e.g. with smoothing turned off) always
%             starts from the same untouched maps.
%
% Name-value options:
%   'SmoothSigmaPx' : sigma (in RVT-map pixels) of a small Gaussian blur applied to each
%                     z's RVT map BEFORE finding its peak (after the RVT transform,
%                     before the trace -- exactly the step you asked for). This exists
%                     because the raw per-pixel argmax can hop between pixels of nearly
%                     equal RVT value frame to frame, adding jitter that has nothing to
%                     do with real particle motion; a very fine blur breaks ties toward
%                     whichever candidate is more robustly a local maximum, without
%                     smearing the map enough to shift a genuine, well-separated peak.
%                     Default 0.75 (px) -- well under the ring-kernel radii (radii_px in
%                     propagate_iscat_rvt_zstack, typically several px), so it denoises
%                     without blurring away real structure. Pass 0 to disable.
%   'Display'     : default false (this is meant to also run unattended over a whole
%                   population, e.g. in a loop building up drift statistics -- you don't
%                   want a figure popping up per particle in that case). If true, opens
%                   a 2x2 figure: x vs. z, y vs. z, the top-down (x,y) path, and a 3D
%                   panel with all three combined -- see "What gets plotted" below.
%   'Axes'        : the 1x4 array of axes handles from a PREVIOUS call's vec.ax, to
%                   overlay a second trace onto the same 4 panels (e.g. sphere vs. bleb).
%                   Default: new figure. The z-vs-position panels are NOT forced to a
%                   common scale with the x-y path panel -- the focus sweep (thousands
%                   of nm) and the lateral drift (often tens of nm) are very different
%                   scales, and forcing them to match would flatten the drift to nothing
%                   (same reasoning applies to the 3D panel's axes).
%   'Color'       : line/marker color, hex string or RGB triplet. Default the dataviz
%                   palette's categorical slot 1 blue (#2a78d6) -- pass slot 2 orange
%                   (#eb6834) for a second overlaid series, this project's usual order.
%   'DisplayName' : legend label. Default: auto-built from the particle's params.
%   'View'        : [az, el] camera view for the 3D panel of a NEW figure, default
%                   [-37.5, 30] (MATLAB's default 3D view).
%
% What gets plotted (each panel uses ONE color per call, so overlaid series stay
% distinguishable -- direction is instead shown by marker shape, not by a color ramp):
%   - a thin line connects consecutive z-planes IN ORDER, with a small filled dot at
%     each z -- this is the actual (smoothed-peak) trace, point by point, not a fit;
%   - the FIRST z is marked with a larger hollow square, the LAST z with a larger
%     filled diamond, so you can read off which way the sweep runs without following
%     the color of a gradient;
%   - any point whose peak pixel falls within one ring radius (rvt_out.radii_px's max)
%     of the image border is drawn HOLLOW instead of filled. Those are flagged because
%     the ring convolution is zero-padded there (see propagate_iscat_rvt_zstack's own
%     border caveat) -- a hollow point sitting out at the edges of a panel is more
%     likely a boundary artifact of the transform than real particle motion, and is
%     worth checking against the raw iSCAT image at that z before reading into it.
%   - the 4th (3D) panel draws the same line/markers in (x,y,z) all at once, so you can
%     see the whole path in one view instead of switching between the three flat panels.
%
% vec fields:
%   x, y, z          : the trace -- x(iz), y(iz) is the (smoothed) RVT-map peak position
%                      (nm) at focus z(iz) (copied straight from rvt_out.z), in order
%   maxRVT           : the (smoothed) RVT value at that peak, at each z (the "center RVT
%                      maximum" you asked to calculate at each focus)
%   peakRow, peakCol : the underlying pixel indices, for reference
%   borderFlag       : true at each z where the peak sits within one ring radius of the
%                      image border (see above) -- these are the hollow points
%   pathLengthNm     : cumulative lateral (x,y) path length across all z steps (sum of
%                      the step-to-step distances) -- this does NOT assume a straight
%                      line; a trace that wanders back and forth has a large
%                      pathLengthNm even though it may end up right back where it started
%   smoothSigmaPx    : the smoothing sigma actually used (0 if disabled)
%   params, mesh_file : copied from rvt_out, for labeling / bookkeeping
%   ax               : 1x4 array of axes handles [x_vs_z, y_vs_z, xy_path, xyz_3d] if
%                      displayed, else []

    p = inputParser;
    addParameter(p, 'SmoothSigmaPx', 0.75);
    addParameter(p, 'Display', false);
    addParameter(p, 'Axes', []);
    addParameter(p, 'Color', []);
    addParameter(p, 'DisplayName', '');
    addParameter(p, 'View', [-37.5, 30]);
    parse(p, varargin{:});
    opt = p.Results;

    assert(isfield(rvt_out, 'rvtMaps') && isfield(rvt_out, 'x'), ...
        'calculate_rvt_vector:noMaps', ...
        ['rvt_out is missing rvtMaps/x -- call propagate_iscat_rvt_zstack with its ' ...
         'default ''StoreMaps'', true (i.e. without ''StoreMaps'', false) first.']);

    x = rvt_out.x;
    z = rvt_out.z(:).';
    nz = size(rvt_out.rvtMaps, 1);
    assert(nz == numel(z), 'calculate_rvt_vector:sizeMismatch', ...
        'rvtMaps has %d z-slices but rvt_out.z has %d elements -- rvt_out looks inconsistent.', ...
        nz, numel(z));

    rmax_px = 1;
    if isfield(rvt_out, 'radii_px') && ~isempty(rvt_out.radii_px)
        rmax_px = max(1, round(max(rvt_out.radii_px)));
    end

    xTraj      = zeros(1, nz);
    yTraj      = zeros(1, nz);
    maxRVT     = zeros(1, nz);
    peakRow    = zeros(1, nz);
    peakCol    = zeros(1, nz);
    borderFlag = false(1, nz);

    for iz = 1:nz
        map = squeeze(rvt_out.rvtMaps(iz, :, :));
        [ny, nx] = size(map);
        if opt.SmoothSigmaPx > 0
            map = gaussian_blur_2d(map, opt.SmoothSigmaPx);   % local copy only -- see docstring
        end
        [v, idx] = max(map(:));
        [r0, c0] = ind2sub(size(map), idx);
        peakRow(iz) = r0;
        peakCol(iz) = c0;
        xTraj(iz) = x(r0);   % dimension 1 = x, matching propagate_iscat_image's convention
        yTraj(iz) = x(c0);   % dimension 2 = y  (same mapping propagate_iscat_rvt_zstack uses)
        maxRVT(iz) = v;
        borderFlag(iz) = (r0 <= rmax_px) || (r0 > ny - rmax_px) || ...
                         (c0 <= rmax_px) || (c0 > nx - rmax_px);
    end

    pathLengthNm = sum(hypot(diff(xTraj), diff(yTraj)));

    vec = struct('x', xTraj, 'y', yTraj, 'z', z, 'maxRVT', maxRVT, ...
        'peakRow', peakRow, 'peakCol', peakCol, 'borderFlag', borderFlag, ...
        'pathLengthNm', pathLengthNm, 'smoothSigmaPx', opt.SmoothSigmaPx, ...
        'params', rvt_out.params, 'mesh_file', rvt_out.mesh_file, 'ax', []);

    % ======================================================================
    % Optional plot: x vs z, y vs z, the top-down (x,y) path, and a 3D panel with all
    % three combined -- the raw (smoothed-peak) trace, nothing fitted or summarized.
    % ======================================================================
    if opt.Display
        blueCol = hex2rgb('#2a78d6');
        gridCol = hex2rgb('#e1e0d9');

        if isempty(opt.Axes)
            fig = figure('Color', 'w', 'Position', [100 100 900 760]);
            tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
            ax = [nexttile(tl), nexttile(tl), nexttile(tl), nexttile(tl)];
            view(ax(4), opt.View(1), opt.View(2));
            for a = ax, hold(a, 'on'); end
        else
            ax = opt.Axes;
            for a = ax, hold(a, 'on'); end
        end

        if isempty(opt.Color)
            col = blueCol;
        elseif ischar(opt.Color) || isstring(opt.Color)
            col = hex2rgb(opt.Color);
        else
            col = opt.Color;
        end
        if isempty(opt.DisplayName)
            name = auto_display_name(rvt_out.params);
        else
            name = opt.DisplayName;
        end

        plot_trace_panel(ax(1), z, xTraj, borderFlag, col, name, 'z focus (nm)', 'x (nm)');
        plot_trace_panel(ax(2), z, yTraj, borderFlag, col, name, 'z focus (nm)', 'y (nm)');
        plot_trace_panel(ax(3), xTraj, yTraj, borderFlag, col, name, 'x (nm)', 'y (nm)');
        plot_trace_panel_3d(ax(4), xTraj, yTraj, z, borderFlag, col, name);
        title(ax(1), 'x vs. focus');
        title(ax(2), 'y vs. focus');
        title(ax(3), 'top-down path');
        title(ax(4), 'all together (3D)');
        sgtitle_txt = 'square = first z, diamond = last z, hollow = near FOV edge';
        try
            title(tl, sgtitle_txt, 'FontWeight', 'normal', 'FontSize', 9, ...
                'Color', [0.34 0.34 0.32]);
        catch
            % no tiledlayout to caption (overlay onto an existing, non-tiled Axes array)
        end

        for i = 1:3
            a = ax(i);
            box(a, 'off');
            a.XColor = [0.34 0.34 0.32];
            a.YColor = [0.34 0.34 0.32];
            a.GridColor = gridCol;
            a.GridAlpha = 1;
            grid(a, 'on');
            a.Layer = 'top';
        end
        a4 = ax(4);
        box(a4, 'off');
        a4.XColor = [0.34 0.34 0.32];
        a4.YColor = [0.34 0.34 0.32];
        a4.ZColor = [0.34 0.34 0.32];
        a4.GridColor = gridCol;
        a4.GridAlpha = 1;
        grid(a4, 'on');
        xlabel(a4, 'x (nm)'); ylabel(a4, 'y (nm)'); zlabel(a4, 'z focus (nm)');

        for a = ax
            nSeries = numel(findobj(a, 'Tag', 'RVTTracePanel'));
            if nSeries > 1
                legend(a, 'show', 'Location', 'best', 'Box', 'off');
            else
                legend(a, 'off');
            end
        end

        vec.ax = ax;
    end
end

function plot_trace_panel(ax, xdat, ydat, borderFlag, col, name, xlab, ylab)
% One 2D panel's worth of "just trace it": a thin connecting line in z_scan order, a
% small filled dot at each point (hollow where borderFlag is true), a larger hollow
% square at the first point and a larger filled diamond at the last.
    plot(ax, xdat, ydat, '-', 'Color', col, 'LineWidth', 1, ...
        'DisplayName', name, 'Tag', 'RVTTracePanel');
    reliable = ~borderFlag;
    if any(reliable)
        plot(ax, xdat(reliable), ydat(reliable), 'o', 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'none', 'MarkerSize', 4, 'HandleVisibility', 'off');
    end
    if any(borderFlag)
        plot(ax, xdat(borderFlag), ydat(borderFlag), 'o', 'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', col, 'MarkerSize', 4, 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    plot(ax, xdat(1), ydat(1), 's', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', col, ...
        'MarkerSize', 9, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(ax, xdat(end), ydat(end), 'd', 'MarkerFaceColor', col, 'MarkerEdgeColor', [1 1 1], ...
        'MarkerSize', 8, 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel(ax, xlab); ylabel(ax, ylab);
end

function plot_trace_panel_3d(ax, xdat, ydat, zdat, borderFlag, col, name)
% The same "just trace it" convention as plot_trace_panel, in 3D (x, y, z together).
    plot3(ax, xdat, ydat, zdat, '-', 'Color', col, 'LineWidth', 1, ...
        'DisplayName', name, 'Tag', 'RVTTracePanel');
    reliable = ~borderFlag;
    if any(reliable)
        plot3(ax, xdat(reliable), ydat(reliable), zdat(reliable), 'o', ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', 'none', 'MarkerSize', 4, ...
            'HandleVisibility', 'off');
    end
    if any(borderFlag)
        plot3(ax, xdat(borderFlag), ydat(borderFlag), zdat(borderFlag), 'o', ...
            'MarkerFaceColor', 'none', 'MarkerEdgeColor', col, 'MarkerSize', 4, ...
            'LineWidth', 1, 'HandleVisibility', 'off');
    end
    plot3(ax, xdat(1), ydat(1), zdat(1), 's', 'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', col, 'MarkerSize', 9, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot3(ax, xdat(end), ydat(end), zdat(end), 'd', 'MarkerFaceColor', col, ...
        'MarkerEdgeColor', [1 1 1], 'MarkerSize', 8, 'LineWidth', 1, 'HandleVisibility', 'off');
end

function out = gaussian_blur_2d(img, sigma_px)
% Separable Gaussian blur, no Image Processing Toolbox dependency (same routine used in
% propagate_iscat_rvt_zstack.m's high-pass step, reused here at a much finer sigma).
    if sigma_px <= 0, out = img; return; end
    r = max(1, ceil(3 * sigma_px));
    t = -r:r;
    g = exp(-(t .^ 2) / (2 * sigma_px ^ 2));
    g = g / sum(g);
    out = conv2(g, g, img, 'same');
end

function name = auto_display_name(params)
    if isfield(params, 'd_core')
        name = sprintf('bleb (core %.0f / bleb %.0f nm)', params.d_core, params.d_bleb);
    else
        lbl = 'sphere';
        if isfield(params, 'label'), lbl = params.label; end
        d = NaN;
        if isfield(params, 'd'), d = params.d; end
        name = sprintf('%s sphere (d = %.0f nm)', lbl, d);
    end
end

function rgb = hex2rgb(h)
    h = char(h);
    if h(1) == '#', h = h(2:end); end
    rgb = [hex2dec(h(1:2)), hex2dec(h(3:4)), hex2dec(h(5:6))] / 255;
end
