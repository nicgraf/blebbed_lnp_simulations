function [outFile, sweep] = render_iscat_zstack_video(particle, sol, specs, z_scan, outFile, varargin)
% RENDER_ISCAT_ZSTACK_VIDEO  Video of the z-by-s iSCAT contrast map as the slicing azimuth
% sweeps through a full rotation, so you can show how asymmetric the iPSF is as the particle
% (or, equivalently, the sampling line through it) rotates.
%
%   render_iscat_zstack_video(particle, sol, specs, z_scan, 'azimuth_zstack.mp4')
%   render_iscat_zstack_video(particle, sol, specs, z_scan, 'out.mp4', ...
%       'AzimuthStep', 10, 'FrameRate', 15)
%
%   particle, sol, specs, z_scan : same as propagate_iscat_zstack / propagate_iscat_azimuth_sweep
%   outFile  : video file path, e.g. 'azimuth_zstack.mp4' (VideoWriter 'MPEG-4' profile)
%
% Name-value options:
%   'AzimuthStep'     : degrees between frames, default 5 (73/72 frames over 0-360)
%   'AzimuthList'     : explicit angle vector, overrides AzimuthStep/IncludeEndpoint
%   'IncludeEndpoint' : include both 0 and 360 (a repeated frame at the loop point,
%                       i.e. a visible pause), default false -- 0:step:355 loops seamlessly
%   'FrameRate'       : video frame rate, default 10
%   'Sweep'           : a precomputed propagate_iscat_azimuth_sweep(...) result, to reuse
%                       across multiple renders (e.g. this video + a report figure) without
%                       recomputing the underlying BEM propagation
%
% The color scale is FIXED across all frames (symmetric about zero, from the sweep's
% global_cmax), not renormalized per frame -- so a frame-to-frame change in contrast is a
% real change in signal, not the color scale adjusting to hide it. A dashed line marks
% best_z (the same z at every azimuth -- see propagate_iscat_azimuth_sweep).
%
% Returns the output file path and the sweep struct (so you can hand `sweep` straight to a
% second render_iscat_zstack_video call, or to render_particle_landing_video's azimuth list,
% without recomputing it).

    p = inputParser;
    addParameter(p, 'AzimuthStep', 5);
    addParameter(p, 'AzimuthList', []);
    addParameter(p, 'IncludeEndpoint', false);
    addParameter(p, 'FrameRate', 10);
    addParameter(p, 'Sweep', []);
    parse(p, varargin{:});

    azList = p.Results.AzimuthList;
    if isempty(azList)
        step = p.Results.AzimuthStep;
        if p.Results.IncludeEndpoint
            azList = 0:step:360;
        else
            azList = 0:step:(360 - step);
        end
    end

    sweep = p.Results.Sweep;
    if isempty(sweep)
        sweep = propagate_iscat_azimuth_sweep(particle, sol, specs, z_scan, azList);
    else
        azList = sweep.azimuth_deg;   % use whatever the precomputed sweep actually covers
    end

    cmax = sweep.global_cmax;
    if cmax == 0, cmax = 1; end
    clims = [-cmax, cmax];

    fig = figure('Color', 'w', 'Position', [100 100 640 640]);
    ax = axes(fig);
    himg = imagesc(ax, sweep.s, sweep.z, sweep.contrast(:, :, 1), clims);
    axis(ax, 'xy');
    colormap(ax, 'gray');
    cb = colorbar(ax);
    cb.Label.String = 'iSCAT contrast (%)';
    xlabel(ax, 's (nm), along azimuth');
    ylabel(ax, 'z focus (nm)');
    hold(ax, 'on');
    yline(ax, sweep.best_z, '--', sprintf('best z = %.0f nm', sweep.best_z), ...
        'Color', [1 0.4 0.4], 'LabelHorizontalAlignment', 'left');
    hold(ax, 'off');
    titleH = title(ax, '');

    vw = VideoWriter(outFile, 'MPEG-4');
    vw.FrameRate = p.Results.FrameRate;
    open(vw);
    cleanupObj = onCleanup(@() close(fig));   % figure closes even if writeVideo errors

    for ia = 1:numel(azList)
        himg.CData = sweep.contrast(:, :, ia);
        titleH.String = sprintf('z-stack, azimuth = %.0f%s  (contrast range \\pm%.2f%%)', ...
            azList(ia), char(176), cmax);
        drawnow;
        writeVideo(vw, getframe(fig));
    end

    close(vw);
    fprintf('wrote %s (%d frames)\n', outFile, numel(azList));
end
