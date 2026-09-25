function out = propagate_iscat_rvt_zstack(mesh_file, specs, z_scan, varargin)
% PROPAGATE_ISCAT_RVT_ZSTACK  Radial-Variance-Transform (RVT) symmetry profile of a
% particle across a z-stack: propagate the iSCAT image at every z, RVT-transform each
% image, and track the value at the most radially-symmetric point (by default the peak
% of the RVT map) as a function of focus height -- so a spherical PSF's symmetry profile
% can be compared directly against a blebbed one, which should be systematically less
% symmetric at defocus because the bleb breaks the particle's rotational symmetry.
%
%   out = propagate_iscat_rvt_zstack(mesh_file, specs, z_scan)
%   out = propagate_iscat_rvt_zstack(mesh_file, specs, z_scan, 'Display', true)
%
%   % overlay a sphere and a bleb on the same axes for direct comparison:
%   outS = propagate_iscat_rvt_zstack('lnp_meshes/full_001_mesh.mat', specs, z_scan);
%   outB = propagate_iscat_rvt_zstack('lnp_meshes/bleb_001_mesh.mat', specs, z_scan, ...
%              'Axes', outS.ax, 'Color', '#eb6834', 'DisplayName', 'bleb 001');
%
%   % also render a video panning the full RVT map through z:
%   out = propagate_iscat_rvt_zstack(mesh_file, specs, z_scan, ...
%              'MakeVideo', true, 'VideoFile', 'bleb_001_rvt.mp4');
%
%   mesh_file : path to a <prefix>_XXX_mesh.mat file that already has a bem_sol appended
%               (by solve_bleb_population / solve_sphere_population). The population type
%               (bleb vs. sphere) is auto-detected from mesh.params (a bleb mesh has a
%               'd_core' field; a sphere mesh's mesh.params.label picks 'full' or 'empty'
%               for load_sphere_solution) -- the right loader is called for you.
%   specs     : make_imaging_specs() output or a prebuilt build_imaging_system() struct
%   z_scan    : focus heights (nm), e.g. -2000:20:2000
%
% ------------------------------------------------------------------------------------
% What "RVT" means here
% ------------------------------------------------------------------------------------
% This is the Radial Variance Transform of Kashkanova, Shkarin, Gholami Mahmoodabadi,
% Blessing, Tuna, Gemeinhardt & Sandoghdar, "Precision single-particle localization
% using radial variance transform," Opt. Express 29, 11070-11083 (2021),
% https://doi.org/10.1364/OE.420670.
%
% The transform itself -- ring-kernel construction and the VoM/MoV variance
% computation -- is NOT reimplemented here. This function calls the project's own
% copy of the original reference MATLAB code (repo folder `rvt/`: RVT.m and
% ringKernel.m, by Anna Kashkanova and Alexey Shkarin, MATLAB translation by Andre
% Gemeinhardt, June 2020) directly, so the numbers this pipeline produces are exactly
% what that verified implementation returns -- no independent reimplementation to
% doubt. This function only does the iSCAT-specific plumbing around it:
%   1) propagate the contrast image at each z;
%   2) optionally high-pass it (img - Gaussian-blur(img, highpass_sigma)) to remove
%      slowly-varying background -- RVT.m itself only mean-subtracts, it does not
%      high-pass, so this step (if enabled) happens before the image is handed to
%      RVT.m;
%   3) optionally upsample the (highpassed) image by an integer factor ('Upsample',
%      default 2) via 2-D spline interpolation, and scale the ring radii by the same
%      factor -- the paper notes this recovers sub-pixel localization precision,
%      because the ring kernels still bin by whole pixels, just finer ones;
%   4) [VoMs, MoVs] = RVT(imgProc, radii_px) -- VoMs is the population variance (across
%      radii) of the per-radius ring-mean images ("basic" RVT / "VoM" -- variance of
%      means -- in the paper); MoVs is the corresponding mean-of-variance term.
%      'Kind' picks which map is actually traced: 'basic' (or its explicit alias
%      'vom') uses VoMs directly -- the MAXIMUM variance of means; 'normalized' uses
%      VoMs./MoVs (the standard RVT ratio, more robust at low SNR / for distorted
%      PSFs, at some cost in subpixel bias -- see the paper); 'mov' skips the RVT
%      ratio entirely and traces the MINIMUM mean of variances (MoVs) instead -- a
%      raw, non-RVT diagnostic, useful for checking whether the VoM/MoV combination
%      itself (rather than the underlying ring convolutions) explains some feature
%      of the trace.
%   5) the peak (maximum for 'basic'/'vom'/'normalized', MINIMUM for 'mov' -- or the
%      geometric-centre point) is located on that fine, upsampled map -- so the
%      trajectory this function tracks (centerRVT, peakXnm, peakYnm) has
%      sub-original-pixel precision -- and only afterward is the map itself resampled
%      back down to the original image resolution (for `rvtMaps`/video, so they stay
%      the same size as `contrastStack` and as before this feature was added).
%      Tracing happens on the fine map; only the stored map is recompressed.
% A radially-symmetric feature (a well-focused, on-axis, non-blebbed PSF) produces a
% sharp, compact peak in the RVT map; breaking that symmetry (defocus asymmetry, a bleb
% off to one side) spreads or lowers the peak. That peak -- 'CenterMode','peak' (default)
% -- is "the most symmetric point" this function tracks across z; 'CenterMode','geometric'
% instead reads the RVT value at the known image centre (sys.icenter), which is only
% meaningful if the particle is centred in the simulated FOV (true for this pipeline).
%
% Requires `rvt/RVT.m` and `rvt/ringKernel.m` (a sibling folder of this file's
% `visualization/` folder, at the repo root) to be on the MATLAB path. This function
% adds it automatically if it can find it relative to its own location; otherwise add
% it yourself, e.g. addpath('<repo root>/rvt').
%
% Name-value options:
%   'RadiusRangeNm'   : [rmin rmax] ring radii in nm, default [100, 1500]. This is the
%                       single most important parameter -- it should bracket the spacing
%                       of the iSCAT interference rings you actually see in the image
%                       (inspect a propagate_iscat_image(...,'Display',true) frame first).
%                       Converted to a plain integer pixel-radius list rmin_px:rmax_px via
%                       sys.pixel_size_nm and passed straight to RVT.m (one thin ring per
%                       integer radius -- no coarsening; see below).
%   'HighpassSigmaNm' : Gaussian high-pass sigma in nm, default 'none' (RVT.m's own
%                       mean-subtraction only, matching a bare call to the reference
%                       code with no extra preprocessing). Pass a number, or 'auto'
%                       (= RadiusRangeNm(2), i.e. rmax) to high-pass before RVT.m runs
%                       if slowly-varying background is swamping the rings. Applied
%                       at the ORIGINAL image resolution, before any upsampling.
%   'Upsample'        : positive integer, default 2. Upsamples the (highpassed) image
%                       by this factor (2-D spline interpolation) before calling
%                       RVT.m, with the ring radii scaled up by the same factor so
%                       they cover the same physical distance -- this is what buys
%                       the sub-pixel precision in centerRVT/peakXnm/peakYnm. Pass 1
%                       to disable (bare original-resolution RVT.m call, no
%                       interpolation anywhere). The stored `rvtMaps` (and the video)
%                       are resampled back down to the original resolution AFTER the
%                       peak is located on the fine map, so they stay the same size
%                       as `contrastStack` regardless of this setting.
%   'Kind'            : which map to trace, and which extremum counts as "the most
%                       symmetric point" (under 'CenterMode','peak') -- see above:
%                         'basic' (default) / 'vom' : VoMs, tracked by its MAXIMUM.
%                         'normalized'               : VoMs./MoVs, tracked by its MAXIMUM.
%                         'mov'                       : MoVs alone (no RVT ratio, and
%                                                        no VoMs computed on top of it),
%                                                        tracked by its MINIMUM.
%   'CenterMode'      : 'peak' (default) or 'geometric' -- see above.
%   'Display'         : default true. Plots centerRVT vs. z on 'Axes' (a new figure if not
%                       given), styled to match this project's other summary figures.
%   'Axes'            : existing axes handle to plot into (for overlaying a second call's
%                       curve on the same axes, e.g. sphere vs. bleb). Default: new figure.
%                       Use the SAME RadiusRangeNm/HighpassSigmaNm/Kind/CenterMode across
%                       overlaid calls -- the curves are only comparable if the transform
%                       was computed the same way for both.
%   'Color'           : line color, hex string or RGB triplet. Default the dataviz palette's
%                       categorical slot 1 blue (#2a78d6) -- pass slot 2 orange (#eb6834)
%                       for a second overlaid series, following this project's fixed
%                       categorical order (see plot_lnp_mesh_stats.m).
%   'DisplayName'     : legend label for this call's curve. Default: auto-built from the
%                       particle's params (e.g. 'bleb (core 70/bleb 50 nm)').
%   'StoreMaps'       : keep the full nz-by-ny-by-nx RVT map stack in `out.rvtMaps` and the
%                       propagated contrast images in `out.contrastStack`, default true.
%                       Forced true if 'MakeVideo' is true (the video needs the maps).
%   'MakeVideo'       : default false. If true, also renders a video panning the full RVT
%                       map through z (fixed color scale, a marker tracking the point used
%                       for centerRVT at each frame).
%   'VideoFile'       : output path, e.g. 'bleb_001_rvt.mp4'. Required if 'MakeVideo' true.
%   'FrameRate'       : video frame rate, default 10.
%
% out fields:
%   z, centerRVT          : the symmetry profile -- the plot's x and y data. centerRVT
%                            is read off the fine, upsampled RVT map (see 'Upsample').
%   peakRow, peakCol       : pixel indices of the point used at each z, INTO THE
%                            UPSAMPLED GRID (size ny*Upsample-by-nx*Upsample; the same
%                            index every z under 'geometric', = the upsampled grid's
%                            centre) -- not indices into rvtMaps/contrastStack, which
%                            stay at the original resolution. Use peakXnm/peakYnm (nm)
%                            for the physical, resolution-independent position.
%   peakXnm, peakYnm       : that point's physical position (nm), per z -- the actual
%                            sub-pixel-precision trajectory this function tracks.
%   rvtMaps, contrastStack : nz-by-ny-by-nx arrays (only if StoreMaps; contrastStack is the
%                            plain iSCAT contrast image at each z, before RVT; rvtMaps is
%                            resampled back down to this same original resolution after
%                            the peak was located on the fine upsampled map)
%   radii_px, radii_nm, highpass_px : the actual ring radii and high-pass sigma used, in
%                            ORIGINAL-resolution pixels/nm (i.e. as applied to rvtMaps'
%                            resolution, not the internal upsampled working resolution)
%   upsample                : the 'Upsample' factor actually used
%   kind, centerMode        : as resolved (after defaulting)
%   mesh_file, particle, sol, params : what was loaded, for reuse without re-loading
%   ax                      : the line-plot axes handle (or [] if Display is false)
%   videoFile               : the video path written (or '' if MakeVideo is false)

    p = inputParser;
    addParameter(p, 'RadiusRangeNm', [100, 1500]);
    addParameter(p, 'HighpassSigmaNm', 'none');
    addParameter(p, 'Upsample', 2);
    addParameter(p, 'Kind', 'basic');
    addParameter(p, 'CenterMode', 'peak');
    addParameter(p, 'Display', true);
    addParameter(p, 'Axes', []);
    addParameter(p, 'Color', []);
    addParameter(p, 'DisplayName', '');
    addParameter(p, 'StoreMaps', true);
    addParameter(p, 'MakeVideo', false);
    addParameter(p, 'VideoFile', '');
    addParameter(p, 'FrameRate', 10);
    parse(p, varargin{:});
    opt = p.Results;

    if ~any(strcmpi(opt.Kind, {'basic', 'vom', 'normalized', 'mov'}))
        error('propagate_iscat_rvt_zstack:badKind', ...
            '''Kind'' must be ''basic'' (or ''vom''), ''normalized'', or ''mov''.');
    end
    useMinExtremum = strcmpi(opt.Kind, 'mov');   % MoV's symmetric point is a MINIMUM,
                                                  % not a maximum like VoM/normalized.
    mapLabel = kind_map_label(lower(opt.Kind));  % for axis/title/colorbar text below
    if ~any(strcmpi(opt.CenterMode, {'peak', 'geometric'}))
        error('propagate_iscat_rvt_zstack:badCenterMode', '''CenterMode'' must be ''peak'' or ''geometric''.');
    end
    if ~isscalar(opt.Upsample) || opt.Upsample < 1 || opt.Upsample ~= round(opt.Upsample)
        error('propagate_iscat_rvt_zstack:badUpsample', ...
            '''Upsample'' must be a positive integer (1 = off).');
    end
    if opt.MakeVideo && isempty(opt.VideoFile)
        error('propagate_iscat_rvt_zstack:noVideoFile', 'Pass ''VideoFile'' when ''MakeVideo'' is true.');
    end
    storeMaps = opt.StoreMaps || opt.MakeVideo;
    if opt.MakeVideo && ~opt.StoreMaps
        warning('propagate_iscat_rvt_zstack:forceStoreMaps', ...
            '''MakeVideo'' needs the RVT maps -- ignoring ''StoreMaps'',false.');
    end

    % ---- make sure the reference RVT.m / ringKernel.m (repo folder `rvt/`) is on the
    % path -- try to find it relative to this file (visualization/../rvt) before failing ----
    if exist('RVT', 'file') ~= 2
        thisFile = mfilename('fullpath');
        rvtDir = fullfile(fileparts(fileparts(thisFile)), 'rvt');
        if exist(fullfile(rvtDir, 'RVT.m'), 'file') == 2
            addpath(rvtDir);
        end
    end
    if exist('RVT', 'file') ~= 2
        error('propagate_iscat_rvt_zstack:noRVT', ...
            ['Can''t find RVT.m on the MATLAB path (expected in <repo root>/rvt, next ' ...
             'to the visualization/ folder). addpath the repo''s rvt/ folder and retry.']);
    end

    % ---- load the particle + solution from the mesh file (population auto-detected) ----
    if ~isfield(specs, 'lens'), specs = build_imaging_system(specs); end
    sys = specs;

    Mpeek = load(mesh_file, 'mesh');
    assert(isfield(Mpeek, 'mesh') && isfield(Mpeek.mesh, 'params'), ...
        'propagate_iscat_rvt_zstack:badFile', ...
        '%s does not look like a mesh file (no mesh.params).', mesh_file);
    isBleb = isfield(Mpeek.mesh.params, 'd_core');
    if isBleb
        [particle, sol] = load_bleb_solution(mesh_file, sys);
    else
        assert(isfield(Mpeek.mesh.params, 'label'), ...
            'propagate_iscat_rvt_zstack:noLabel', ...
            '%s has no mesh.params.label to pick ''full''/''empty'' for load_sphere_solution.', ...
            mesh_file);
        [particle, sol] = load_sphere_solution(mesh_file, sys, Mpeek.mesh.params.label);
    end

    far = farfields(sol, sys.lens.dir);   % once per particle, focus-independent
    x = sys.x;
    px_nm = sys.pixel_size_nm;
    nz = numel(z_scan);

    % ---- image size from a throwaway first frame ----
    contrast0 = propagate_iscat_image(sol, z_scan(1), sys, far);
    [ny, nx] = size(contrast0);

    % ---- ring radii (pixels) -- a plain integer list, passed straight to RVT.m -- and
    % the high-pass sigma ----
    rmin_px = max(1, round(opt.RadiusRangeNm(1) / px_nm));
    rmax_px = max(rmin_px + 1, round(opt.RadiusRangeNm(2) / px_nm));
    if 2 * rmax_px + 1 > min(ny, nx)
        warning('propagate_iscat_rvt_zstack:largeRadius', ...
            ['RadiusRangeNm(2) (%.0f nm = %d px) makes the ring kernel (%d px across) as ' ...
             'large as the image (%d x %d px) -- results near rmax will be dominated by ' ...
             'zero-padding. Consider a smaller RadiusRangeNm(2).'], ...
            opt.RadiusRangeNm(2), rmax_px, 2 * rmax_px + 1, ny, nx);
    end
    radii_px = rmin_px:rmax_px;
    radii_nm = radii_px * px_nm;

    % ---- upsampled working grid (for sub-pixel peak localization -- see 'Upsample'
    % in the docstring). Ring radii are scaled by the same integer factor so they
    % still cover the same physical distance on the finer grid. ----
    upsample = round(opt.Upsample);
    doUpsample = upsample > 1;
    if doUpsample
        nx_up = nx * upsample;
        x_up = linspace(x(1), x(end), nx_up);
        radii_px_up = radii_px * upsample;
        [~, icenter_up] = min(abs(x_up));
    else
        nx_up = nx;
        x_up = x;
        radii_px_up = radii_px;
        icenter_up = sys.icenter;
    end

    if ischar(opt.HighpassSigmaNm) || isstring(opt.HighpassSigmaNm)
        switch lower(char(opt.HighpassSigmaNm))
            case 'auto', hp_px = rmax_px;
            case 'none', hp_px = [];
            otherwise
                error('propagate_iscat_rvt_zstack:badHighpass', ...
                    '''HighpassSigmaNm'' must be ''auto'', ''none'', or a number (nm).');
        end
    elseif isempty(opt.HighpassSigmaNm)
        hp_px = [];
    else
        hp_px = max(0, opt.HighpassSigmaNm) / px_nm;
    end

    % ---- main loop: propagate -> (optional highpass) -> (optional upsample) ->
    % RVT.m -> locate the trajectory point on the FINE map -> recompress the map ----
    centerRVT = zeros(nz, 1);
    peakRow   = zeros(nz, 1);
    peakCol   = zeros(nz, 1);
    peakXnm   = zeros(nz, 1);
    peakYnm   = zeros(nz, 1);
    if storeMaps
        rvtMaps       = zeros(nz, ny, nx);
        contrastStack = zeros(nz, ny, nx);
    end

    for iz = 1:nz
        if iz == 1
            img = contrast0;
        else
            img = propagate_iscat_image(sol, z_scan(iz), sys, far);
        end
        if storeMaps, contrastStack(iz, :, :) = img; end

        imgProc = img;
        if ~isempty(hp_px) && hp_px > 0
            imgProc = imgProc - gaussian_blur_2d(imgProc, hp_px);
        end
        if doUpsample
            imgForRVT = resample_image_2d(imgProc, x, nx_up, 'spline');
        else
            imgForRVT = imgProc;
        end
        % RVT.m mean-subtracts internally -- hand it the (optionally high-passed,
        % optionally upsampled) image directly, no separate mean-subtraction here.
        [VoMs, MoVs] = RVT(imgForRVT, radii_px_up);

        switch lower(opt.Kind)
            case {'basic', 'vom'}
                rvtMapFine = VoMs;
            case 'normalized'
                rvtMapFine = VoMs ./ MoVs;
            case 'mov'
                rvtMapFine = MoVs;
        end

        % Trace on the fine (upsampled) map -- this is what gives centerRVT and the
        % peak position sub-original-pixel precision. 'mov' tracks the MINIMUM
        % (useMinExtremum); every other Kind tracks the MAXIMUM.
        switch lower(opt.CenterMode)
            case 'peak'
                if useMinExtremum
                    [v, idx] = min(rvtMapFine(:));
                else
                    [v, idx] = max(rvtMapFine(:));
                end
                [r0, c0] = ind2sub(size(rvtMapFine), idx);
            case 'geometric'
                r0 = icenter_up; c0 = icenter_up;
                v = rvtMapFine(r0, c0);
        end
        centerRVT(iz) = v;
        peakRow(iz) = r0;
        peakCol(iz) = c0;
        peakXnm(iz) = x_up(r0);
        peakYnm(iz) = x_up(c0);

        % Only now recompress the map back down to the original resolution, so
        % rvtMaps/video stay the same size as contrastStack regardless of 'Upsample'.
        if storeMaps
            if doUpsample
                rvtMaps(iz, :, :) = resample_image_2d(rvtMapFine, x_up, nx, 'spline');
            else
                rvtMaps(iz, :, :) = rvtMapFine;
            end
        end
    end

    % ---- output struct ----
    out = struct('z', z_scan(:).', 'centerRVT', centerRVT(:).', ...
        'peakRow', peakRow(:).', 'peakCol', peakCol(:).', ...
        'peakXnm', peakXnm(:).', 'peakYnm', peakYnm(:).', ...
        'x', x, ...   % pixel coordinate grid (nm), shared by both dims -- lets a
                       % downstream function (e.g. calculate_rvt_vector) turn any
                       % rvtMaps pixel index back into a physical position without
                       % needing to rebuild sys/specs.
        'radii_px', radii_px, 'radii_nm', radii_nm, 'highpass_px', hp_px, ...
        'upsample', upsample, ...
        'kind', lower(opt.Kind), 'centerMode', lower(opt.CenterMode), ...
        'mesh_file', mesh_file, 'particle', particle, 'sol', sol, 'params', particle.params, ...
        'ax', [], 'videoFile', '');
    if storeMaps
        out.rvtMaps = rvtMaps;
        out.contrastStack = contrastStack;
    end

    % ======================================================================
    % Display: centerRVT vs z, styled to match this project's other figures
    % ======================================================================
    if opt.Display
        col.blue = hex2rgb('#2a78d6');
        gridCol  = hex2rgb('#e1e0d9');

        if isempty(opt.Axes)
            fig = figure('Color', 'w', 'Position', [100 100 640 460]);
            ax = axes(fig);
        else
            ax = opt.Axes;
        end
        hold(ax, 'on');

        if isempty(opt.Color)
            lineCol = col.blue;
        elseif ischar(opt.Color) || isstring(opt.Color)
            lineCol = hex2rgb(opt.Color);
        else
            lineCol = opt.Color;
        end
        if isempty(opt.DisplayName)
            name = auto_display_name(particle.params, isBleb);
        else
            name = opt.DisplayName;
        end

        plot(ax, z_scan, centerRVT, '-', 'Color', lineCol, 'LineWidth', 2, ...
            'DisplayName', name, 'Tag', 'RVTSeries');
        if useMinExtremum
            [~, ib] = min(centerRVT);
        else
            [~, ib] = max(centerRVT);
        end
        plot(ax, z_scan(ib), centerRVT(ib), 'o', 'MarkerFaceColor', lineCol, ...
            'MarkerEdgeColor', [1 1 1], 'MarkerSize', 7, 'HandleVisibility', 'off');

        style_axes(ax, gridCol);
        xlabel(ax, 'z focus (nm)');
        if strcmpi(opt.CenterMode, 'geometric')
            pointDesc = 'image-centre';
        elseif useMinExtremum
            pointDesc = 'min-peak';
        else
            pointDesc = 'max-peak';
        end
        ylabel(ax, sprintf('%s at %s point (a.u.)', mapLabel, pointDesc));
        title(ax, sprintf('%s vs. focus', mapLabel));

        nSeries = numel(findobj(ax, 'Tag', 'RVTSeries'));
        if nSeries > 1
            legend(ax, 'show', 'Location', 'best', 'Box', 'off');
        else
            legend(ax, 'off');
        end

        out.ax = ax;
    end

    % ======================================================================
    % Optional video: RVT map panning through z
    % ======================================================================
    if opt.MakeVideo
        cmax = max(rvtMaps(:));
        if cmax <= 0, cmax = 1; end
        clims = [0, cmax];
        seqCmap = seq_blue_colormap(256);

        vfig = figure('Color', 'w', 'Position', [100 100 640 620]);
        vax = axes(vfig);
        himg = imagesc(vax, x, x, squeeze(rvtMaps(1, :, :)), clims);
        axis(vax, 'image');
        colormap(vax, seqCmap);
        cb = colorbar(vax);
        cb.Label.String = sprintf('%s (a.u.)', mapLabel);
        xlabel(vax, 'x (nm)'); ylabel(vax, 'y (nm)');
        hold(vax, 'on');
        markH = plot(vax, peakXnm(1), peakYnm(1), 'o', 'MarkerEdgeColor', [1 1 1], ...
            'MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1.5);
        hold(vax, 'off');
        titleH = title(vax, '');

        vw = VideoWriter(opt.VideoFile, 'MPEG-4');
        vw.FrameRate = opt.FrameRate;
        open(vw);
        cleanupObj = onCleanup(@() close(vfig));

        for iz = 1:nz
            himg.CData = squeeze(rvtMaps(iz, :, :));
            markH.XData = peakXnm(iz);
            markH.YData = peakYnm(iz);
            titleH.String = sprintf('%s map, z = %.0f nm  (%s point marked)', ...
                mapLabel, z_scan(iz), lower(opt.CenterMode));
            drawnow;
            writeVideo(vw, getframe(vfig));
        end
        close(vw);
        fprintf('wrote %s (%d frames)\n', opt.VideoFile, nz);
        out.videoFile = opt.VideoFile;
    end
end

% ==========================================================================
% Local helpers
% ==========================================================================

function imgOut = resample_image_2d(imgIn, xIn, nOut, method)
% Resample a square image imgIn (on a symmetric grid xIn, shared by both dims) onto a
% new nOut-by-nOut grid spanning the SAME physical extent [xIn(1), xIn(end)], via 2-D
% interpolation. Used both ways: to upsample the image before RVT.m (nOut > numel(xIn),
% for sub-pixel peak localization) and to recompress the resulting RVT map back down to
% the original resolution afterward (nOut < numel(xIn)) -- same routine either way, so
% the two operations are exact inverses of each other's grid geometry.
    xOut = linspace(xIn(1), xIn(end), nOut);
    [Xin, Yin] = meshgrid(xIn, xIn);
    [Xout, Yout] = meshgrid(xOut, xOut);
    imgOut = interp2(Xin, Yin, imgIn, Xout, Yout, method);
end

function out = gaussian_blur_2d(img, sigma_px)
% Separable Gaussian blur, no Image Processing Toolbox dependency. Used only for the
% optional high-pass preprocessing step before RVT.m -- the RVT computation itself
% (ring kernels, VoM/MoV) is entirely delegated to rvt/RVT.m + rvt/ringKernel.m.
    if sigma_px <= 0, out = img; return; end
    r = max(1, ceil(3 * sigma_px));
    t = -r:r;
    g = exp(-(t .^ 2) / (2 * sigma_px ^ 2));
    g = g / sum(g);
    out = conv2(g, g, img, 'same');
end

function name = auto_display_name(params, isBleb)
    if isBleb
        name = sprintf('bleb (core %.0f / bleb %.0f nm)', params.d_core, params.d_bleb);
    else
        lbl = 'sphere';
        if isfield(params, 'label'), lbl = params.label; end
        name = sprintf('%s sphere (d = %.0f nm)', lbl, params.d);
    end
end

function lbl = kind_map_label(kindLower)
% Human-readable label for the map/quantity being traced, used in axis labels, plot
% titles, and the video colorbar -- keeps that text accurate for 'mov' (which is
% deliberately NOT "RVT") as well as the RVT-proper kinds.
    switch kindLower
        case {'basic', 'vom'}
            lbl = 'VoM (variance of means)';
        case 'normalized'
            lbl = 'Normalized RVT (VoM / MoV)';
        case 'mov'
            lbl = 'MoV (mean of variances)';
        otherwise
            lbl = kindLower;
    end
end

function style_axes(ax, gridCol)
    box(ax, 'off');
    ax.XColor = [0.34 0.34 0.32];
    ax.YColor = [0.34 0.34 0.32];
    ax.GridColor = gridCol;
    ax.GridAlpha = 1;
    grid(ax, 'on');
    ax.Layer = 'top';
end

function cmap = seq_blue_colormap(n)
% Dataviz-skill sequential "blue" ramp (light -> dark, steps 100-700), interpolated to
% n rows, for the RVT map (a non-negative, magnitude-only quantity -- sequential, not
% diverging or rainbow, per the skill's color-formula rules).
    if nargin < 1, n = 256; end
    hexSteps = {'#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', ...
                '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b'};
    rgb = cellfun(@hex2rgb, hexSteps, 'UniformOutput', false);
    rgb = cat(1, rgb{:});
    xOld = linspace(0, 1, size(rgb, 1));
    xNew = linspace(0, 1, n);
    cmap = interp1(xOld, rgb, xNew, 'pchip');
    cmap = min(max(cmap, 0), 1);
end

function rgb = hex2rgb(h)
    h = char(h);
    if h(1) == '#', h = h(2:end); end
    rgb = [hex2dec(h(1:2)), hex2dec(h(3:4)), hex2dec(h(5:6))] / 255;
end
