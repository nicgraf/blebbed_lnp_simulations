function sweep = propagate_iscat_azimuth_sweep(particle, sol, specs, z_scan, azimuth_deg_list)
% PROPAGATE_ISCAT_AZIMUTH_SWEEP  z-by-s iSCAT stack at MANY slicing azimuths, computed
% efficiently: the underlying x-y contrast image at each z does not depend on azimuth (only
% which line through it gets sampled does), so this builds that image ONCE per z and slices
% it at every requested angle, instead of calling propagate_iscat_zstack once per angle
% (which would redundantly recompute the same x-y images n_azimuth times over).
%
%   sweep = propagate_iscat_azimuth_sweep(particle, sol, specs, z_scan, 0:5:355)
%
%   particle          : struct from load_bleb_solution or load_sphere_solution
%   sol               : solved stratified.solution for that particle
%   specs             : make_imaging_specs() output or a prebuilt build_imaging_system() struct
%   z_scan            : focus heights (nm)
%   azimuth_deg_list  : vector of slicing angles (deg, from +x toward +y)
%
% sweep fields:
%   contrast, numerator, ibg  : nz-by-ns-by-n_azimuth arrays
%   s (nm), z (nm), azimuth_deg (the input vector, as given)
%   best_z   : SCALAR, not per-azimuth -- the s=0 sample is the same physical point (the
%              particle centre) at every azimuth, so the focus-selection rule (same one
%              propagate_iscat_zstack uses: peak |numerator| at centre, excluding z where
%              the reference intensity has collapsed below specs.ibg_floor_frac of its
%              scan max) gives an identical answer regardless of slicing angle.
%   global_cmax : max(abs(contrast(:))) over the WHOLE sweep, for a consistent color scale
%                 across azimuths (e.g. for a video) -- unlike propagate_iscat_zstack's
%                 per-call cmax, which would rescale every frame and hide the asymmetry.
%   params   : copied from particle

    if ~isfield(specs, 'lens'), specs = build_imaging_system(specs); end
    sys = specs;
    x = sys.x;  s = x;
    az = azimuth_deg_list(:).';        % row
    naz = numel(az);
    nz = numel(z_scan);  ns = numel(s);

    % all azimuths' query points at once, per z (ns-by-naz)
    Xq = s(:) * cosd(az);
    Yq = s(:) * sind(az);

    far = farfields(sol, sys.lens.dir);      % once for the whole sweep

    contrast  = zeros(nz, ns, naz);
    numerator = zeros(nz, ns, naz);
    ibgArr    = zeros(nz, ns, naz);

    for iz = 1:nz
        [~, im, ibg] = propagate_iscat_image(sol, z_scan(iz), sys, far);
        if ~sys.dim1_is_x, im = im.'; ibg = ibg.'; end   % make dim 1 = x
        num = im - ibg;

        Fn = griddedInterpolant({x, x}, num, 'linear', 'none');
        Fb = griddedInterpolant({x, x}, ibg, 'linear', 'none');
        Nq = Fn(Xq, Yq);   % ns-by-naz, one interpolant build+eval per z for ALL azimuths
        Bq = Fb(Xq, Yq);

        numerator(iz, :, :) = reshape(Nq, [1, ns, naz]);
        ibgArr(iz, :, :)    = reshape(Bq, [1, ns, naz]);
        contrast(iz, :, :)  = reshape(Nq ./ Bq * 100, [1, ns, naz]);
    end

    % focus selection at s = 0 -- identical physical point regardless of azimuth, so this
    % is computed once (any azimuth slice gives the same numbers, up to interpolation noise)
    [~, ic] = min(abs(s));
    num0 = squeeze(numerator(:, ic, 1));
    ibg0 = squeeze(ibgArr(:, ic, 1));
    ok = ibg0 >= sys.ibg_floor_frac * max(ibg0);
    score = abs(num0);
    if all(~ok)
        warning('propagate_iscat_azimuth_sweep:allExcluded', ...
            'Every z has ibg below the floor -- widen z_scan or lower ibg_floor_frac.');
    else
        score(~ok) = -Inf;
    end
    [~, ib] = max(score);
    best_z = z_scan(ib);

    sweep = struct('contrast', contrast, 'numerator', numerator, 'ibg', ibgArr, ...
                   's', s, 'z', z_scan(:).', 'azimuth_deg', az, 'best_z', best_z, ...
                   'global_cmax', max(abs(contrast(:))), 'params', particle.params);
end
