function specs = make_imaging_specs()
% MAKE_IMAGING_SPECS  Single source of truth for the optics + camera used for BOTH the
% BEM solve and the iSCAT propagation, so the two can't drift out of sync.
% Edit values here; everything else takes `specs` (or the `sys` built from it).

    % ---- source ----
    specs.lambda = 555;            % nm  (update if your acquisition channel differs)
    specs.pol    = [1, 0, 0];      % single linear polarization (see note in the readme text)
    specs.dir    = [0, 0, 1];      % epi illumination: from the glass side, upward

    % ---- media ----
    specs.n_glass   = 1.515;
    specs.n_water   = 1.335;       % also the medium around the particles
    specs.n_lipid   = 1.495;       % Leslie & Cullis
    specs.n_aqueous = 1.45;        % placeholder
    specs.n_full    = 1.515;       % full (mRNA-loaded) LNP sphere population
    specs.n_empty   = 1.335;       % empty LNP sphere population; placeholder -- same as
                                   % n_water for now, i.e. no core contrast, update if a
                                   % better estimate exists
    specs.idx_water = 2; specs.idx_lipid = 3; specs.idx_aqueous = 4;   % positions in mat_set
    specs.idx_full  = 5; specs.idx_empty = 6;                          % (spheres are appended,
                                                                       % not renumbered, so
                                                                       % existing bleb solves
                                                                       % stay valid)

    % ---- sphere mesh (full/empty populations) ----
    specs.nverts = 144;   % trisphere discretization -- bump for larger particles / higher
                          % index contrast; refine until cross sections stop changing

    % ---- substrate gap (shared with the bleb mesh builder) ----
    specs.gap = 2;   % nm, height of each particle's LOWEST point above z=0

    % ---- BEM ----
    specs.bem_order = 5;
    specs.tau_method = 'native';   % 'native' | 'merge' | 'concat'(debug only); see rebuild_bleb_tau

    % ---- objective + magnification (Nikon Ti2, 100x plan apo NA 1.45, 1.5x intermediate) ----
    specs.NA               = 1.45;
    specs.mag_objective    = 100;
    specs.mag_intermediate = 1.5;
    specs.mag_total        = specs.mag_objective * specs.mag_intermediate;   % informational

    % ---- camera (Prime95B). 118 um FOV over 1608 px already includes the 150x system
    %      magnification, so this is the effective pixel pitch AT THE SAMPLE. ----
    specs.sensor_fov_um = 118;
    specs.sensor_pixels = 1608;
    specs.roi_pixels_full = 1200;  % documented experimental ROI, for reference only
    specs.npix = 61;               % simulated window per particle (odd, centered on the particle)

    % ---- image array convention ----
    % efield(...) is assumed to return arrays whose dimension 1 is x and dimension 2 is y
    % (as implied by the original flip(...,1) "flip x-axis" step). Set false if your images
    % turn out transposed.
    specs.dim1_is_x = true;

    % ---- z-focus selection (used only to report best_z in the z-stack) ----
    specs.ibg_floor_frac = 0.1;    % ignore z where reference intensity < 10% of its scan max
end
