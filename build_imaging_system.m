function sys = build_imaging_system(specs)
% BUILD_IMAGING_SYSTEM  Turn the plain `specs` struct into everything the solver and the
% imaging functions need: materials, layer, incident field, lens, reference field, pixel grid.
% Build this ONCE and reuse it (propagate_* accept either specs or sys).

    sys = specs;
    sys.k0 = 2 * pi / specs.lambda;

    sys.mat_glass   = Material(specs.n_glass^2,   1);
    sys.mat_water   = Material(specs.n_water^2,   1);
    sys.mat_lipid   = Material(specs.n_lipid^2,   1);
    sys.mat_aqueous = Material(specs.n_aqueous^2, 1);
    sys.mat_full    = Material(specs.n_full^2,    1);
    sys.mat_empty   = Material(specs.n_empty^2,   1);
    % layer materials first (glass, water), per the nanobem convention -- every population
    % (bleb, full, empty) shares this ONE mat_set so material indices are never ambiguous;
    % full/empty are appended at the end so existing bleb solves (indices 1-4) stay valid.
    sys.mat_set = [sys.mat_glass, sys.mat_water, sys.mat_lipid, sys.mat_aqueous, ...
                   sys.mat_full, sys.mat_empty];

    sys.layer = stratified.layerstructure([sys.mat_glass, sys.mat_water], 0);
    sys.einc  = optics.decompose(sys.k0, specs.pol, specs.dir);

    sys.pixel_size_nm = specs.sensor_fov_um * 1000 / specs.sensor_pixels;
    n = specs.npix;
    sys.x = sys.pixel_size_nm * (-(n-1)/2 : (n-1)/2);   % nm, centered on the particle
    sys.icenter = ceil(n / 2);

    % imaging lens: object side = glass, image side = air, epi detection (rot = 180 about y)
    air = Material(1, 1);
    rot = optics.roty(180);
    sys.lens = optics.lensimage(sys.mat_glass, air, sys.k0, specs.NA, 'rot', rot);

    % iSCAT reference wave = reflection of the excitation at the substrate, heading back down
    einc = sys.einc;
    sys.refl = secondary(einc, sys.layer, 'dir', 'down');
end
