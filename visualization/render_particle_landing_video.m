function outFile = render_particle_landing_video(particle, specs, outFile, varargin)
% RENDER_PARTICLE_LANDING_VIDEO  Companion video to render_iscat_zstack_video: the particle
% mesh sitting on the glass, held FIXED, with a line through its centre that sweeps through
% the same azimuth angles -- so a viewer can see exactly what "rotating the azimuth" means
% physically, alongside the corresponding z-stack video.
%
%   render_particle_landing_video(particle, specs, 'azimuth_landing.mp4')
%   render_particle_landing_video(particle, specs, 'out.mp4', 'AzimuthStep', 10)
%   render_particle_landing_video(particle, specs, 'out.mp4', 'AzimuthList', sweep.azimuth_deg)
%
%   particle : struct from load_bleb_solution or load_sphere_solution (.params, .tau).
%              Population type is auto-detected the same way as propagate_iscat_zstack:
%              a bleb particle has particle.params.d_core, a sphere has particle.params.d.
%   specs    : make_imaging_specs() output OR a prebuilt build_imaging_system() struct --
%              either works, only the idx_lipid/idx_aqueous/idx_water/idx_full/idx_empty
%              fields are used (no BEM/optics setup needed for a geometry-only video).
%   outFile  : video file path, e.g. 'azimuth_landing.mp4'
%
% Name-value options:
%   'AzimuthStep'     : degrees between frames, default 5
%   'AzimuthList'      : explicit angle vector -- pass sweep.azimuth_deg from
%                        render_iscat_zstack_video/propagate_iscat_azimuth_sweep to make the
%                        two videos frame-matched
%   'IncludeEndpoint'  : include both 0 and 360, default false (seamless loop)
%   'FrameRate'        : video frame rate, default 10
%   'View'             : [az, el] camera view, default [40, 22] (MATLAB's `view` convention)
%   'Orbit'             : if true, the CAMERA also slowly spins once around over the course
%                         of the video (purely cosmetic, off by default so the only thing
%                         that visibly moves is the azimuth line, matching what this video
%                         is meant to show)
%
% For a bleb: lipid core (gold) and aqueous bleb (blue) are colored separately and a black
% arrow marks the fixed core->bleb axis (from the actual solved mesh, not recomputed from
% angles). For a sphere: the whole mesh in one color (blue = 'full', orange = 'empty', by
% whichever material it was actually solved against; gray otherwise). Both show the glass
% plane at z=0 and a red dashed line through the particle centre that sweeps through the
% requested azimuths, with a marker at its +s end (the direction render_iscat_zstack_video's
% s > 0 points toward).

    p = inputParser;
    addParameter(p, 'AzimuthStep', 5);
    addParameter(p, 'AzimuthList', []);
    addParameter(p, 'IncludeEndpoint', false);
    addParameter(p, 'FrameRate', 10);
    addParameter(p, 'View', [40, 22]);
    addParameter(p, 'Orbit', false);
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

    isBleb = isfield(particle.params, 'd_core');
    tau = particle.tau;
    pos = vertcat(tau.pos);
    xyRange = max(pos(:, 1:2), [], 1) - min(pos(:, 1:2), [], 1);
    pad = 0.25 * max(xyRange);
    xlims = [min(pos(:,1)) - pad, max(pos(:,1)) + pad];
    ylims = [min(pos(:,2)) - pad, max(pos(:,2)) + pad];

    fig = figure('Color', 'w', 'Position', [100 100 720 640]);
    ax = axes(fig);
    hold(ax, 'on');

    % glass surface at z = 0
    patch(ax, xlims([1 2 2 1]), ylims([1 1 2 2]), [0 0 0 0], [0.75 0.85 0.95], ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none');

    if isBleb
        io = vertcat(tau.inout);
        maskLipidOut   = io(:, 1) == specs.idx_lipid   & io(:, 2) == specs.idx_water;
        maskAqueousOut = io(:, 1) == specs.idx_aqueous & io(:, 2) == specs.idx_water;
        maskInterface  = ~maskLipidOut & ~maskAqueousOut;

        if any(maskLipidOut)
            plot(tau(maskLipidOut),   'FaceColor', [0.95 0.70 0.20], 'FaceAlpha', 0.95);
        end
        if any(maskAqueousOut)
            plot(tau(maskAqueousOut), 'FaceColor', [0.30 0.60 0.90], 'FaceAlpha', 0.95);
        end
        if any(maskInterface)
            plot(tau(maskInterface),  'FaceColor', [0.60 0.60 0.60], 'FaceAlpha', 0.4);
        end

        if any(maskLipidOut) && any(maskAqueousOut)
            coreCenter = mean(pos(maskLipidOut, :), 1);
            blebCenter = mean(pos(maskAqueousOut, :), 1);
            quiver3(ax, coreCenter(1), coreCenter(2), coreCenter(3), ...
                blebCenter(1) - coreCenter(1), blebCenter(2) - coreCenter(2), ...
                blebCenter(3) - coreCenter(3), 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.8);
        else
            coreCenter = mean(pos, 1);
        end

        titleBase = sprintf('tilt = %.0f%s, gold = lipid core, blue = aqueous bleb', ...
            particle.params.theta_tilt_deg, char(176));
    else
        io = tau(1).inout;
        matName = 'particle';  matColor = [0.55 0.55 0.55];
        if isfield(specs, 'idx_full') && io(1) == specs.idx_full
            matName = 'full';   matColor = [0.16 0.47 0.84];
        elseif isfield(specs, 'idx_empty') && io(1) == specs.idx_empty
            matName = 'empty';  matColor = [0.92 0.41 0.20];
        end
        plot(tau, 'FaceColor', matColor, 'FaceAlpha', 0.95);
        coreCenter = mean(pos, 1);
        titleBase = sprintf('%s sphere, d = %.1f nm', matName, particle.params.d);
    end

    % rotating slicing-axis line (updated per frame) + a marker at its +s end
    L = 0.6 * max(xyRange);
    lineH = plot3(ax, nan(1,2), nan(1,2), nan(1,2), '--', 'Color', [0.85 0.2 0.2], 'LineWidth', 2);
    markH = plot3(ax, nan, nan, nan, 'o', 'MarkerFaceColor', [0.85 0.2 0.2], ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 8);

    hold(ax, 'off');
    axis(ax, 'equal');
    grid(ax, 'on');
    view(ax, p.Results.View(1), p.Results.View(2));
    xlabel(ax, 'x (nm)'); ylabel(ax, 'y (nm)'); zlabel(ax, 'z (nm)');
    titleH = title(ax, '');

    vw = VideoWriter(outFile, 'MPEG-4');
    vw.FrameRate = p.Results.FrameRate;
    open(vw);
    cleanupObj = onCleanup(@() close(fig));

    n = numel(azList);
    for ia = 1:n
        az = azList(ia);
        dirv = [cosd(az), sind(az), 0];
        lineH.XData = coreCenter(1) + [-L, L] * dirv(1);
        lineH.YData = coreCenter(2) + [-L, L] * dirv(2);
        lineH.ZData = coreCenter(3) * [1, 1];
        markH.XData = coreCenter(1) + L * dirv(1);
        markH.YData = coreCenter(2) + L * dirv(2);
        markH.ZData = coreCenter(3);

        if p.Results.Orbit
            view(ax, p.Results.View(1) + 360 * (ia - 1) / n, p.Results.View(2));
        end

        titleH.String = sprintf('%s  |  slicing azimuth = %.0f%s', titleBase, az, char(176));
        drawnow;
        writeVideo(vw, getframe(fig));
    end

    close(vw);
    fprintf('wrote %s (%d frames)\n', outFile, n);
end
