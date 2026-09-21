function out = propagate_iscat_zstack(particle, sol, azimuth_deg, specs, z_scan)
% PROPAGATE_ISCAT_ZSTACK  z-by-s iSCAT stack along the in-plane axis at angle azimuth_deg.
%
%   out = propagate_iscat_zstack(particle, sol, azimuth_deg, specs, z_scan)
%
%   particle    : struct from load_bleb_solution (only .params is used, for metadata)
%   sol         : solved stratified.solution for that particle
%   azimuth_deg : in-plane angle of the slicing axis, measured from +x toward +y. This is
%                 the neck-axis azimuth used when the mesh was built (particle.params.
%                 phi_azimuth_deg); pass [] to use that value.
%   specs       : make_imaging_specs() output or a prebuilt build_imaging_system() struct
%   z_scan      : focus heights (nm), e.g. -2000:20:2000
%
% The x-y image is computed at every z and sampled along the line through the particle
% centre at the given azimuth, so s > 0 points toward the bleb side and s < 0 away from it.
% The far field is computed once, since it does not depend on focus.
%
% out fields (nz-by-ns arrays have rows = z, columns = s):
%   contrast, numerator (im - ibg, raw interference), ibg, s (nm), z (nm),
%   azimuth_deg, best_z, params (copied from particle)
%
% best_z: focus where |numerator| at s = 0 peaks, ignoring z where the reference intensity
% has collapsed below specs.ibg_floor_frac of its scan maximum (same rule as before).

    if isempty(azimuth_deg)
        azimuth_deg = particle.params.phi_azimuth_deg;
    end
    if ~isfield(specs, 'lens'), specs = build_imaging_system(specs); end
    sys = specs;

    x = sys.x;
    s = x;                                      % same pitch and extent as the image axes
    xq = s * cosd(azimuth_deg);
    yq = s * sind(azimuth_deg);

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
        numerator(iz, :) = Fn(xq, yq);
        ibg_line(iz, :)  = Fb(xq, yq);
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

    out = struct('contrast', contrast, 'numerator', numerator, 'ibg', ibg_line, ...
                 's', s, 'z', z_scan(:).', 'azimuth_deg', azimuth_deg, ...
                 'best_z', z_scan(ib), 'params', particle.params);
end