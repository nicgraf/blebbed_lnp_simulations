function out = propagate_iscat_zstack(particle, sol, specs, z_scan, varargin)
% PROPAGATE_ISCAT_ZSTACK  z-by-s iSCAT stack along the in-plane axis at angle azimuth_deg.
%
%   out = propagate_iscat_zstack(particle, sol, specs, z_scan)
%   out = propagate_iscat_zstack(particle, sol, specs, z_scan, 'Display', true)
%   out = propagate_iscat_zstack(particle, sol, specs, z_scan, 'AzimuthDeg', 30, 'Display', true)
%
%   particle : struct from load_bleb_solution OR load_sphere_solution (.params for
%              metadata, .tau for the mesh geometry, used only for the landing-view
%              display). Population type is auto-detected from particle.params: a bleb
%              particle has a 'd_core' field, a sphere particle has a 'd' field.
%   sol      : solved stratified.solution for that particle
%   specs    : make_imaging_specs() output or a prebuilt build_imaging_system() struct
%   z_scan   : focus heights (nm), e.g. -2000:20:2000
%
% Name-value options:
%   'AzimuthDeg' : in-plane angle of the slicing axis, measured from +x toward +y.
%                  OPTIONAL -- if you don't pass it, it's picked automatically:
%                    bleb   -> -particle.params.phi_azimuth_deg (the mesh's own neck azimuth)
%                    sphere -> 0 (arbitrary; a sphere has no natural azimuth, it's
%                               rotationally symmetric, so any angle gives the same stack)
%                  Pass a number to override either case (e.g. to compare a bleb and a
%                  sphere sliced at the same nominal angle).
%   'Display'    : default false. If true, opens TWO figures:
%                    1) the z-by-s contrast map, colormap scaled symmetrically about zero
%                       (gray = no signal), with a marker at best_z.
%                    2) a 3D rendering of the particle mesh as it landed on the glass.
%                       For a bleb: lipid core and aqueous bleb colored separately, plus an
%                       arrow along the core->bleb axis. For a sphere: the whole mesh in
%                       one color (blue for a 'full' population, orange for 'empty', by
%                       whichever material index it was actually solved against -- gray if
%                       neither matches). Both also show the glass plane at z=0 and the
%                       in-plane slicing axis used for the stack above.
%
% The x-y image is computed at every z and sampled along the line through the particle
% centre at the given azimuth, so s > 0 points toward the bleb side and s < 0 away from it
% (for a sphere, s > 0/< 0 are equivalent by symmetry).
% The far field is computed once, since it does not depend on focus.
%
% out fields (nz-by-ns arrays have rows = z, columns = s):
%   contrast, numerator (im - ibg, raw interference), ibg, s (nm), z (nm),
%   azimuth_deg, best_z, params (copied from particle),
%   ax (z-stack axes handle if displayed, else []), mesh_ax (landing-view axes, else [])
%
% best_z: focus where |numerator| at s = 0 peaks, ignoring z where the reference intensity
% has collapsed below specs.ibg_floor_frac of its scan maximum (same rule as before).

    isBleb = isfield(particle.params, 'd_core');

    p = inputParser;
    if isBleb
        defaultAzimuth = particle.params.phi_azimuth_deg;
    else
        defaultAzimuth = 0;   % arbitrary for a sphere -- no natural azimuth to fall back to
    end
    addParameter(p, 'AzimuthDeg', defaultAzimuth);
    addParameter(p, 'Display', false);
    parse(p, varargin{:});
    azimuth_deg = p.Results.AzimuthDeg;
    if isempty(azimuth_deg), azimuth_deg = defaultAzimuth; end   % tolerate azimuth_deg=[]
    displayIt = p.Results.Display;

    if ~isfield(specs, 'lens'), specs = build_imaging_system(specs); end
    sys = specs;
    x = sys.x;
    s = x;                                      % same pitch and extent as the image axes
    xq = (s * cosd(azimuth_deg)).';
    yq = (s * sind(azimuth_deg)).';
    far = farfields(sol, sys.lens.dir);         % once per particle
    nz = numel(z_scan);  ns = numel(s);
    contrast  = zeros(nz, ns);
    numerator = zeros(nz, ns);
    ibg_line  = zeros(nz, ns);
    for iz = 1:nz
        [~, im, ibg] = propagate_iscat_image(sol, z_scan(iz), sys, far);
        num = im - ibg;
        if ~sys.dim1_is_x, num = num.'; ibg = ibg.'; end   % make dim 1 = x
        Fn = griddedInterpolant({x, x}, num, 'linear', 'none');
        Fb = griddedInterpolant({x, x}, ibg, 'linear', 'none');
        numerator(iz, :) = Fn(xq, yq).';
        ibg_line(iz, :)  = Fb(xq, yq).';
    end
    contrast = numerator ./ ibg_line * 100;

    % focus selection at the particle centre (s = 0), excluding collapsed-reference z
    [~, ic] = min(abs(s));
    ibg_c = ibg_line(:, ic);
    ok = ibg_c >= sys.ibg_floor_frac * max(ibg_c);
    score = abs(numerator(:, ic));
    if all(~ok)
        warning('propagate_iscat_zstack:allExcluded', ...
            'Every z has ibg below the floor -- widen z_scan or lower ibg_floor_frac.');
    else
        score(~ok) = -Inf;
    end
    [~, ib] = max(score);
    best_z = z_scan(ib);

    ax = [];
    mesh_ax = [];
    if displayIt
        % ---- figure 1: z-by-s contrast map ----
        cmax = max(abs(contrast(:)));
        if cmax == 0, cmax = 1; end   % avoid a degenerate [0 0] clim on an all-zero stack
        clims = [-cmax, cmax];

        figure;
        ax = axes;
        imagesc(ax, s, z_scan, contrast, clims);
        axis(ax, 'xy');               % z increases upward
        colormap(ax, 'gray');         % symmetric gray LUT: 0 contrast = mid-gray
        cb = colorbar(ax);
        cb.Label.String = 'iSCAT contrast (%)';
        xlabel(ax, 's (nm), along azimuth');
        ylabel(ax, 'z focus (nm)');

        hold(ax, 'on');
        yline(ax, best_z, '--', sprintf('best z = %.0f nm', best_z), ...
            'Color', [1 0.4 0.4], 'LabelHorizontalAlignment', 'left');
        hold(ax, 'off');

        title(ax, sprintf('%s, azimuth = %.0f%s  (contrast range \\pm%.2f%%)', ...
            'z-stack', azimuth_deg, char(176), cmax));

        % ---- figure 2: how the particle landed on the glass ----
        mesh_ax = render_particle_landing(particle, sys, azimuth_deg, isBleb);
    end

    out = struct('contrast', contrast, 'numerator', numerator, 'ibg', ibg_line, ...
                 's', s, 'z', z_scan(:).', 'azimuth_deg', azimuth_deg, ...
                 'best_z', best_z, 'params', particle.params, 'ax', ax, 'mesh_ax', mesh_ax);
end

function ax2 = render_particle_landing(particle, sys, azimuth_deg, isBleb)
% RENDER_PARTICLE_LANDING  3D view of the solved boundary mesh sitting on the glass (z=0).
% Bleb: lipid core and aqueous bleb colored separately, so you can see the landing
% tilt/azimuth and the shape directly -- not a re-derivation from params, so it can't drift
% out of sync with what was actually meshed and solved.
% Sphere: the whole (homogeneous) mesh in one color, chosen from whichever material index
% it was actually solved against (sys.idx_full / sys.idx_empty), not from params.label, so
% this stays correct even if a mesh's label metadata is missing or stale.

    tau = particle.tau;
    pos = vertcat(tau.pos);
    xyRange = max(pos(:, 1:2), [], 1) - min(pos(:, 1:2), [], 1);
    pad = 0.25 * max(xyRange);
    xlims = [min(pos(:,1)) - pad, max(pos(:,1)) + pad];
    ylims = [min(pos(:,2)) - pad, max(pos(:,2)) + pad];

    figure;
    ax2 = axes;
    hold(ax2, 'on');

    % glass surface at z = 0
    patch(ax2, xlims([1 2 2 1]), ylims([1 1 2 2]), [0 0 0 0], [0.75 0.85 0.95], ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none');

    if isBleb
        io = vertcat(tau.inout);   % Nx2, [insideMatIdx, outsideMatIdx] per element
        maskLipidOut   = io(:, 1) == sys.idx_lipid   & io(:, 2) == sys.idx_water;
        maskAqueousOut = io(:, 1) == sys.idx_aqueous & io(:, 2) == sys.idx_water;
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

        % core->bleb axis, straight from the mesh (robust to any tilt/azimuth convention,
        % since it is measured from the actual solved geometry, not recomputed from angles)
        if any(maskLipidOut) && any(maskAqueousOut)
            coreCenter = mean(pos(maskLipidOut, :), 1);
            blebCenter = mean(pos(maskAqueousOut, :), 1);
            quiver3(ax2, coreCenter(1), coreCenter(2), coreCenter(3), ...
                blebCenter(1) - coreCenter(1), blebCenter(2) - coreCenter(2), ...
                blebCenter(3) - coreCenter(3), 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.8);
        else
            coreCenter = mean(pos, 1);
        end

        titleStr = sprintf(['landing view: tilt = %.0f%s, azimuth = %.0f%s  ' ...
            '(gold = lipid core, blue = aqueous bleb, gray plane = glass)'], ...
            particle.params.theta_tilt_deg, char(176), azimuth_deg, char(176));
    else
        % single homogeneous boundary -- classify by the material it was actually solved
        % against, not by params.label, so this can't silently mislabel a stale/missing tag
        io = tau(1).inout;
        matName = 'particle';  matColor = [0.55 0.55 0.55];
        if isfield(sys, 'idx_full') && io(1) == sys.idx_full
            matName = 'full';   matColor = [0.16 0.47 0.84];   % dataviz blue
        elseif isfield(sys, 'idx_empty') && io(1) == sys.idx_empty
            matName = 'empty';  matColor = [0.92 0.41 0.20];   % dataviz orange
        end
        plot(tau, 'FaceColor', matColor, 'FaceAlpha', 0.95);

        coreCenter = mean(pos, 1);
        d = particle.params.d;
        titleStr = sprintf('landing view: %s sphere, d = %.1f nm  (gray plane = glass)', ...
            matName, d);
    end

    % in-plane slicing axis used for the z-stack above (the s-axis), drawn through the
    % particle centroid at its own height for reference
    L = 0.6 * max(xyRange);
    dir = [cosd(azimuth_deg), sind(azimuth_deg), 0];
    plot3(ax2, coreCenter(1) + [-L, L] * dir(1), coreCenter(2) + [-L, L] * dir(2), ...
        coreCenter(3) * [1 1], '--', 'Color', [0.85 0.2 0.2], 'LineWidth', 1.5);

    hold(ax2, 'off');
    axis(ax2, 'equal');
    grid(ax2, 'on');
    view(ax2, 40, 22);
    xlabel(ax2, 'x (nm)'); ylabel(ax2, 'y (nm)'); zlabel(ax2, 'z (nm)');
    title(ax2, titleStr);
end
