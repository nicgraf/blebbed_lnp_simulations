function tau = rebuild_bleb_tau(mesh, mat_set, idx_water, idx_lipid, idx_aqueous, method)
% REBUILD_BLEB_TAU  Turn a lean mesh struct (from build_fused_bleb_mesh) back into the
% boundary `tau` that the BEM solver / stratified.solution need.
%
% method:
%   'native' (default) : ONE call BoundaryEdge(mat, p1, inout1, p2, inout2, p3, inout3),
%                        i.e. nanobem's own multi-patch constructor, which numbers the
%                        edge basis functions globally so the patches are properly linked.
%   'merge'            : build the three patches separately, then re-number through
%                        ShapeEdge (what merge_bleb_patches in ISCATSimulatorApp does).
%   'concat'           : plain concatenation. NOT valid for solving (each patch restarts
%                        its edge numbering at 1); kept only for debugging.
%
% Solving and loading MUST use the same method (it fixes the element/edge order that the
% saved e and h refer to). tau is returned as a ROW array, which is what nanobem's
% own code (e.g. BoundaryElement/touching) expects.

    if nargin < 6 || isempty(method), method = 'native'; end

    node = mesh.node;
    p_lipid_out   = particle(node, double(mesh.f_lipid_out));
    p_aqueous_out = particle(node, double(mesh.f_aqueous_out));
    p_interface   = particle(node, double(mesh.f_interface));

    io_lipid_out   = [idx_lipid,   idx_water];     % lipid   | water
    io_aqueous_out = [idx_aqueous, idx_water];     % aqueous | water
    io_interface   = [idx_lipid,   idx_aqueous];   % lipid   | aqueous

    switch method
        case 'native'
            tau = BoundaryEdge(mat_set, p_lipid_out, io_lipid_out, ...
                                        p_aqueous_out, io_aqueous_out, ...
                                        p_interface, io_interface);
        case {'merge', 'concat'}
            t1 = BoundaryEdge(mat_set, p_lipid_out,   io_lipid_out);
            t2 = BoundaryEdge(mat_set, p_aqueous_out, io_aqueous_out);
            t3 = BoundaryEdge(mat_set, p_interface,   io_interface);
            tau_all = [t1(:).', t2(:).', t3(:).'];
            if strcmp(method, 'concat')
                tau = tau_all;
            else
                tau = BoundaryEdge(ShapeEdge(tau_all));
                [~, i1] = ismember(vertcat(tau_all.pos), vertcat(tau.pos), 'rows');
                tau = tau(i1);
            end
        otherwise
            error('rebuild_bleb_tau:method', 'Unknown method ''%s''.', method);
    end
    tau = tau(:).';
end
