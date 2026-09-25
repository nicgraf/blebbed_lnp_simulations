function [contrast, im, ibg, ax] = propagate_iscat_image(sol, z, specs, far, displayIt)
% PROPAGATE_ISCAT_IMAGE  x-by-y iSCAT image of a solved particle at focus height z (nm).
%
%   [contrast, im, ibg] = propagate_iscat_image(sol, z, specs)
%   [contrast, im, ibg] = propagate_iscat_image(sol, z, sys, far)          % reuse far fields
%   [contrast, im, ibg, ax] = propagate_iscat_image(sol, z, specs, [], true)  % also display
%
%   sol       : stratified.solution (e.g. from load_bleb_solution)
%   z         : focus position (nm); a scalar, same convention as focus = [0 0 z] before
%   specs     : camera/imaging specs from make_imaging_specs, or the prebuilt system from
%               build_imaging_system (recommended inside loops -- avoids rebuilding the lens)
%   far       : optional, precomputed farfields(sol, lens.dir) -- it does not depend on z.
%               Pass [] to let this function compute it.
%   displayIt : optional, default false. If true, opens a figure showing the contrast
%               image with a colormap scaled symmetrically about zero (gray = no signal),
%               since iSCAT contrast is signed (constructive vs. destructive interference).
%
%   contrast : iSCAT contrast in percent, (im - ibg)./ibg*100, npix-by-npix
%   im       : total interference intensity |E_sca + E_ref|^2
%   ibg      : reference-only intensity |E_ref|^2
%   ax       : axes handle of the displayed image (only created if displayIt is true;
%              [] otherwise)
%
% Axes: dimension 1 = x, dimension 2 = y, at sys.x (nm, pixel pitch at the sample).

    if nargin < 4, far = []; end
    if nargin < 5 || isempty(displayIt), displayIt = false; end

    if ~isfield(specs, 'lens'), specs = build_imaging_system(specs); end
    sys = specs;
    lens = sys.lens; refl = sys.refl; x = sys.x;

    if isempty(far)
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

    ax = [];
    if displayIt
        cmax = max(abs(contrast(:)));
        if cmax == 0, cmax = 1; end   % avoid a degenerate [0 0] clim on an all-zero image
        clims = [-cmax, cmax];

        figure;
        ax = axes;
        imagesc(ax, x, x, contrast, clims);
        axis(ax, 'image');
        colormap(ax, 'gray');         % symmetric gray LUT: 0 contrast = mid-gray
        cb = colorbar(ax);
        cb.Label.String = 'iSCAT contrast (%)';
        xlabel(ax, 'x (nm)'); ylabel(ax, 'y (nm)');
        title(ax, sprintf('z = %.0f nm  (contrast range \\pm%.2f%%)', z, cmax));
    end
end