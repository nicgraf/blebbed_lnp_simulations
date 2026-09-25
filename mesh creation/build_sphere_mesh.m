function mesh = build_sphere_mesh(d, nverts, gap)
% BUILD_SPHERE_MESH  Lean mesh struct for a single homogeneous sphere (full or empty LNP),
% same storage philosophy as build_fused_bleb_mesh.m: geometry + particle parameters only,
% no material, no BoundaryEdge/particle objects (those get rebuilt cheaply at solve/load
% time via rebuild_sphere_tau.m).
%
%   d      : sphere diameter (nm)
%   nverts : trisphere discretization (see make_imaging_specs -- refine until cross
%            sections stop changing)
%   gap    : height (nm) of the sphere's LOWEST point above the substrate interface z=0
%
% mesh.params : d, nverts, gap
% mesh.node   : Nx3 double, vertices (already shifted so the sphere sits `gap` above z=0)
% mesh.faces  : Mx3 int32

    p = trisphere(nverts, d);
    p.verts(:, 3) = p.verts(:, 3) + 0.5 * d + gap;   % lowest point sits `gap` above z=0

    mesh.node   = p.verts;
    mesh.faces  = int32(p.faces);
    mesh.params = struct('d', d, 'nverts', nverts, 'gap', gap);
end
