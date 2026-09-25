function mesh = build_fused_bleb_mesh(d_core, d_bleb, neck_overlap, step, radbound, gap, ...
        theta_tilt_deg, phi_azimuth_deg)
% BUILD_FUSED_BLEB_MESH  Lean version. Returns a small struct instead of nanobem objects.
%
% Extracted verbatim from lnp_mesh_creation.mlx into its own file so it can be
% reused (e.g. by build_bleb_population_grid.m) without copy-pasting it out of the
% live script each time. The body below is UNCHANGED from the live script --
% geometry, orientation convention, tet handling and ELEMENT ORDER are identical to
% the original (so solutions saved from the old meshes still line up).
%
% The struct holds only what is needed to rebuild the boundary (see
% rebuild_bleb_tau.m) plus the particle parameters, so per-bleb files are tiny
% and can be told apart later:
%
%   mesh.params  : d_core, d_bleb, neck_overlap, theta_tilt_deg, phi_azimuth_deg,
%                  gap, step, radbound
%   mesh.node    : Nx3 double, surface vertices only, already rotated/shifted
%   mesh.f_lipid_out, mesh.f_aqueous_out, mesh.f_interface : Mx3 int32 faces
%
% Materials (mat_set, idx_*) are not needed here; pass them to rebuild_bleb_tau.
%
% Requires iso2mesh on the path (for vol2mesh).

    r_core = d_core / 2;
    r_bleb = d_bleb / 2;
    center_dist = r_core + r_bleb - neck_overlap;
    pad = 3 * step;
    x = (-r_core - pad):step:(center_dist + r_bleb + pad);
    y = (-max(r_core, r_bleb) - pad):step:(max(r_core, r_bleb) + pad);
    z = y;
    [X, Y, Z] = ndgrid(x, y, z);
    in_core = (X.^2 + Y.^2 + Z.^2) <= r_core^2;
    in_bleb = ((X - center_dist).^2 + Y.^2 + Z.^2) <= r_bleb^2;
    img = zeros(size(X), 'uint8');
    img(in_core) = 1;               % lipid
    img(in_bleb & ~in_core) = 2;    % aqueous (lipid wins the overlap)
    clear X Y Z in_core in_bleb

    opt.radbound = radbound;
    [node, elem] = vol2mesh(img, 1:size(img,1), 1:size(img,2), 1:size(img,3), ...
        opt, 1000, 1, 'cgalmesh');
    clear img
    node = node(:, 1:3);

    % rotate neck axis to (theta_tilt, phi_azimuth), then sit `gap` above z = 0
    theta = deg2rad(theta_tilt_deg);
    phi   = deg2rad(phi_azimuth_deg);
    Ry90  = [0 0 1; 0 1 0; -1 0 0];
    Rtilt = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    Razim = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    node = (Razim * Rtilt * Ry90 * node')';
    node(:, 3) = node(:, 3) - min(node(:, 3)) + gap;

    tets   = elem(:, 1:4);
    labels = elem(:, 5);
    nt = size(tets, 1);
    if numel(unique(labels)) < 2
        error('build_fused_bleb_mesh:notFused', ...
            ['Only one labeled region was meshed -- the aqueous bleb was likely fully ' ...
             'engulfed by the lipid core or too thin to resolve at this step. Increase ' ...
             'd_bleb relative to neck_overlap, decrease neck_overlap, or use a finer step.']);
    end

    % consistent (positive-volume) tet orientation
    v1 = node(tets(:,1), :); v2 = node(tets(:,2), :);
    v3 = node(tets(:,3), :); v4 = node(tets(:,4), :);
    flip = dot(cross(v2 - v1, v3 - v1, 2), v4 - v1, 2) < 0;
    tmp = tets(flip, 3); tets(flip, 3) = tets(flip, 4); tets(flip, 4) = tmp;
    clear v1 v2 v3 v4

    % outward-oriented faces of each tet
    localFaces = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
    allFaces  = zeros(nt * 4, 3);
    faceLabel = zeros(nt * 4, 1);
    for i = 1:4
        allFaces((i-1)*nt+1:i*nt, :) = tets(:, localFaces(i, :));
        faceLabel((i-1)*nt+1:i*nt)   = labels;
    end

    % vectorised face classification (same result AND order as the old per-group loop)
    [~, ~, ic] = unique(sort(allFaces, 2), 'rows');
    n      = accumarray(ic, 1);            n = n(ic);       % tets sharing this face
    labsum = accumarray(ic, faceLabel);    labsum = labsum(ic);
    pick = @(mask) allFaces(sortByGroup(find(mask), ic), :);
    f_lo = pick(n == 1 & faceLabel == 1);                       % lipid | medium
    f_ao = pick(n == 1 & faceLabel == 2);                       % aqueous | medium
    f_if = pick(n == 2 & labsum == 3 & faceLabel == 1);         % lipid | aqueous, lipid-side orientation

    % keep only vertices actually used by surface faces
    used = unique([f_lo(:); f_ao(:); f_if(:)]);
    map = zeros(size(node, 1), 1);
    map(used) = 1:numel(used);
    mesh.node = node(used, :);

    mesh.f_lipid_out   = int32(map(f_lo));
    mesh.f_aqueous_out = int32(map(f_ao));
    mesh.f_interface   = int32(map(f_if));

    mesh.params = struct('d_core', d_core, 'd_bleb', d_bleb, 'neck_overlap', neck_overlap, ...
        'theta_tilt_deg', theta_tilt_deg, 'phi_azimuth_deg', phi_azimuth_deg, ...
        'gap', gap, 'step', step, 'radbound', radbound);

    fprintf('  fused bleb mesh: %d + %d + %d = %d elements\n', size(f_lo,1), size(f_ao,1), ...
        size(f_if,1), size(f_lo,1) + size(f_ao,1) + size(f_if,1));
end

function rows = sortByGroup(rows, ic)
    [~, o] = sort(ic(rows));
    rows = rows(o);
end
