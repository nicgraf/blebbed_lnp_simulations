function mesh = adjust_bleb_azimuth(mesh, new_phi_azimuth_deg)
% ADJUST_BLEB_AZIMUTH  Change a saved bleb mesh's azimuth angle WITHOUT rebuilding the
% mesh. build_fused_bleb_mesh (voxelize + vol2mesh) is the slow step -- this is a
% single 3x3 rotation applied to the stored vertices, effectively instant.
%
%   mesh2 = adjust_bleb_azimuth(mesh, new_phi_azimuth_deg)
%
%   mesh                 : a mesh struct from build_fused_bleb_mesh /
%                          build_bleb_population_grid (needs mesh.node and
%                          mesh.params.phi_azimuth_deg)
%   new_phi_azimuth_deg  : the azimuth angle (degrees) to rotate to
%
% Why this is exact, not an approximation: in build_fused_bleb_mesh, the core
% sphere is centered at the ORIGIN in x,y for every theta_tilt/phi_azimuth (all
% three rotations -- Ry90, Rtilt, Razim -- are about the origin, so they can't move
% it off x=y=0); only the LAST step, a shift of node(:,3) by `gap - min(z)`, moves
% the particle up to sit `gap` above the substrate, and a rotation about z leaves
% every z-coordinate (and therefore that shift) unchanged. Razim(phi_azimuth) -- a
% plain rotation about z -- is applied to the whole node set BEFORE that shift.
% So re-azimuthing to a new angle is exactly an ADDITIONAL rotation about z by
% (new - old) degrees: it reproduces, to floating-point precision, what calling
% build_fused_bleb_mesh again with the new phi_azimuth_deg (everything else
% unchanged) would give -- no re-voxelization, re-meshing, or re-shift needed.
% mesh.f_lipid_out / f_aqueous_out / f_interface (face topology) don't change at
% all, since rotation doesn't touch which vertices make up which face.
%
% THETA_TILT_DEG is NOT this simple: the z-shift itself depends on min(z), which
% DOES change under a change of tilt (a different tilt changes which point of the
% fused shape is lowest), so an existing mesh's tilt can't be adjusted this way --
% that needs a full rebuild via build_fused_bleb_mesh with the new theta_tilt_deg.

    if ~isfield(mesh, 'node') || ~isfield(mesh, 'params') || ...
            ~isfield(mesh.params, 'phi_azimuth_deg')
        error('adjust_bleb_azimuth:badMesh', ...
            'mesh must have mesh.node and mesh.params.phi_azimuth_deg (a build_fused_bleb_mesh output).');
    end

    delta = deg2rad(new_phi_azimuth_deg - mesh.params.phi_azimuth_deg);
    Rz = [cos(delta) -sin(delta) 0; sin(delta) cos(delta) 0; 0 0 1];
    mesh.node = (Rz * mesh.node')';
    mesh.params.phi_azimuth_deg = mod(new_phi_azimuth_deg, 360);
end
