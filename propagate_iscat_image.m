function [contrast, im, ibg] = propagate_iscat_image(sol, z, specs, far)
% PROPAGATE_ISCAT_IMAGE  x-by-y iSCAT image of a solved particle at focus height z (nm).
%
%   [contrast, im, ibg] = propagate_iscat_image(sol, z, specs)
%   [contrast, im, ibg] = propagate_iscat_image(sol, z, sys, far)   % reuse far fields
%
%   sol   : stratified.solution (e.g. from load_bleb_solution)
%   z     : focus position (nm); a scalar, same convention as focus = [0 0 z] before
%   specs : camera/imaging specs from make_imaging_specs, or the prebuilt system from
%           build_imaging_system (recommended inside loops -- avoids rebuilding the lens)
%   far   : optional, precomputed farfields(sol, lens.dir) -- it does not depend on z
%
%   contrast : iSCAT contrast in percent, (im - ibg)./ibg*100, npix-by-npix
%   im       : total interference intensity |E_sca + E_ref|^2
%   ibg      : reference-only intensity |E_ref|^2
%
% Axes: dimension 1 = x, dimension 2 = y, at sys.x (nm, pixel pitch at the sample).

    if ~isfield(specs, 'lens'), specs = build_imaging_system(specs); end
    sys = specs;
    lens = sys.lens; refl = sys.refl; x = sys.x;

    if nargin < 4 || isempty(far)
        far = farfields(sol, lens.dir);
    end

    focus = [0, 0, z];
    isca = efield(lens, far,  x, x, 'focus', focus);
    iref = efield(lens, refl, x, x, 'focus', focus);
    isca = flip(isca, 1);            % consequence of the rot = 180 lens axis
    iref = flip(iref, 1);

    im       = dot(isca + iref, isca + iref, 3);
    ibg      = dot(iref, iref, 3);
    contrast = (im - ibg) ./ ibg * 100;
end