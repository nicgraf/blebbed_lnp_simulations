function stats = plot_lnp_mesh_stats(sphere_mesh_dir, bleb_mesh_dir, varargin)
% PLOT_LNP_MESH_STATS  Summary figure of spherical- and blebbed-LNP mesh statistics.
%
%   stats = plot_lnp_mesh_stats(sphere_mesh_dir, bleb_mesh_dir)
%   stats = plot_lnp_mesh_stats({full_dir, empty_dir}, bleb_mesh_dir)   % multiple sphere folders
%   stats = plot_lnp_mesh_stats(sphere_mesh_dir, bleb_mesh_dir, 'NumBins', 25, ...
%                                'SaveAs', 'lnp_mesh_stats.png')
%
%   sphere_mesh_dir : folder (or cell array of folders) holding <prefix>_XXX_mesh.mat
%                     files from build_sphere_population.m (e.g. 'full_001_mesh.mat',
%                     'empty_001_mesh.mat'). Only mesh.params is read -- fast even if
%                     bem_sol has already been appended.
%   bleb_mesh_dir   : folder (or cell array) holding bleb_XXX_mesh.mat files from the
%                     blebbed-LNP population builder.
%
% Name-value options:
%   'SpherePattern' : glob for sphere files, default '*_mesh.mat'
%   'BlebPattern'   : glob for bleb files,   default 'bleb_*_mesh.mat'
%   'NumBins'       : histogram bin count, default 30
%   'SaveAs'        : path to save the figure (png/pdf/...); [] (default) = don't save
%
% Geometry used for the blebbed-LNP derived metrics (matches build_fused_bleb_mesh.m):
% the two fused spheres sit center_dist = r_core + r_bleb - neck_overlap apart, so the
% overall particle length along the neck axis is
%     total_length = r_core + center_dist + r_bleb = d_core + d_bleb - neck_overlap
% and the aspect ratio is that length divided by the larger of the two sphere diameters
% (the widest cross-section either lobe reaches):
%     aspect_ratio = total_length / max(d_core, d_bleb)
%
% Output `stats` is a struct with the raw per-particle vectors (sphere_diameter_nm by
% label, d_core_nm, d_bleb_nm, neck_overlap_nm, total_length_nm, aspect_ratio, plus the
% bleb orientation vectors theta_tilt_deg / phi_azimuth_deg) plus mean/std/n for each, so
% the numbers behind the figure are available without re-parsing.

    p = inputParser;
    addParameter(p, 'SpherePattern', '*_mesh.mat');
    addParameter(p, 'BlebPattern', 'bleb_*_mesh.mat');
    addParameter(p, 'NumBins', 30);
    addParameter(p, 'SaveAs', []);
    parse(p, varargin{:});
    nbins = p.Results.NumBins;

    % ---- categorical palette (dataviz skill default; fixed order, not cycled) ----
    col.blue   = hex2rgb('#2a78d6');
    col.orange = hex2rgb('#eb6834');
    col.aqua   = hex2rgb('#1baf7a');
    col.yellow = hex2rgb('#eda100');
    edgeCol    = [1 1 1];
    gridCol    = hex2rgb('#e1e0d9');

    % ---- gather sphere params (grouped by mesh.params.label, e.g. 'full'/'empty') ----
    sphereParams = collect_mesh_params(sphere_mesh_dir, p.Results.SpherePattern);
    if isempty(sphereParams)
        warning('plot_lnp_mesh_stats:noSpheres', 'No sphere mesh files found in %s.', ...
            strjoin(cellstr_of(sphere_mesh_dir), ', '));
    end
    if ~isempty(sphereParams) && isfield(sphereParams, 'label')
        sphereLabels = {sphereParams.label};
    else
        sphereLabels = repmat({'sphere'}, size(sphereParams));
    end
    groupNames = unique(sphereLabels, 'stable');   % encounter order = fixed color order

    % ---- gather bleb params ----
    blebParams = collect_mesh_params(bleb_mesh_dir, p.Results.BlebPattern);
    if isempty(blebParams)
        warning('plot_lnp_mesh_stats:noBlebs', 'No bleb mesh files found in %s.', ...
            strjoin(cellstr_of(bleb_mesh_dir), ', '));
    end

    d_core = [blebParams.d_core];
    d_bleb = [blebParams.d_bleb];
    neck_overlap = [blebParams.neck_overlap];
    total_length = d_core + d_bleb - neck_overlap;
    largest_width = max(d_core, d_bleb);
    aspect_ratio = total_length ./ largest_width;

    % Orientation. By construction (isotropic sampling: theta_tilt_deg = acosd(2*rand-1),
    % phi_azimuth_deg = 360*rand):
    %   tilt is a POLAR angle (0-180 deg) -- NOT uniform in angle, it is uniform over the
    %     sphere, so its density is sin(theta)/2 (theta in radians); plain linear mean/std
    %     are still valid summary stats since 0/180 deg are true endpoints, not a wrap point.
    %   azimuth is a CIRCULAR angle (0-360 deg, wraps) -- naive linear mean/std are misleading
    %     near the wrap, so circular statistics are used instead.
    theta_tilt_deg = [blebParams.theta_tilt_deg];
    phi_azimuth_deg = [blebParams.phi_azimuth_deg];

    C = mean(cosd(phi_azimuth_deg));
    S = mean(sind(phi_azimuth_deg));
    circMeanAzimuth = mod(atan2d(S, C), 360);
    circR = sqrt(C.^2 + S.^2);
    circStdAzimuth = sqrt(-2 * log(circR)) * 180 / pi;

    % ======================================================================
    % Figure
    % ======================================================================
    fig = figure('Color', 'w', 'Position', [100 100 1560 720]);
    tl = tiledlayout(fig, 2, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, 'LNP mesh population statistics', 'FontWeight', 'bold', 'FontSize', 13);

    % ---- (1) spherical LNP diameters, grouped by label ----
    ax1 = nexttile(tl);
    groupCols = {col.blue, col.orange, col.aqua, col.yellow};
    hold(ax1, 'on');
    legEntries = gobjects(0);
    for g = 1:numel(groupNames)
        d_g = [sphereParams(strcmp(sphereLabels, groupNames{g})).d];
        h = histogram(ax1, d_g, nbins, 'FaceColor', groupCols{mod(g-1,4)+1}, ...
            'FaceAlpha', 0.85, 'EdgeColor', edgeCol, 'LineWidth', 0.5);
        legEntries(end+1) = h; %#ok<AGROW>
    end
    style_axes(ax1, gridCol);
    xlabel(ax1, 'diameter (nm)'); ylabel(ax1, 'count');
    title(ax1, sprintf('Spherical LNP diameters (n = %d)', numel(sphereLabels)));
    if numel(groupNames) > 1
        legend(ax1, legEntries, groupNames, 'Location', 'best', 'Box', 'off');
    end
    hold(ax1, 'off');

    % ---- (2) core diameter ----
    ax2 = nexttile(tl);
    histogram(ax2, d_core, nbins, 'FaceColor', col.blue, 'FaceAlpha', 0.85, ...
        'EdgeColor', edgeCol, 'LineWidth', 0.5);
    style_axes(ax2, gridCol);
    xlabel(ax2, 'core (lipid) diameter (nm)'); ylabel(ax2, 'count');
    title(ax2, sprintf('Bleb core diameter (n = %d, %.1f\\pm%.1f nm)', ...
        numel(d_core), mean(d_core), std(d_core)));

    % ---- (3) bleb diameter ----
    ax3 = nexttile(tl);
    histogram(ax3, d_bleb, nbins, 'FaceColor', col.orange, 'FaceAlpha', 0.85, ...
        'EdgeColor', edgeCol, 'LineWidth', 0.5);
    style_axes(ax3, gridCol);
    xlabel(ax3, 'bleb (aqueous) diameter (nm)'); ylabel(ax3, 'count');
    title(ax3, sprintf('Bleb diameter (n = %d, %.1f\\pm%.1f nm)', ...
        numel(d_bleb), mean(d_bleb), std(d_bleb)));

    % ---- (4) overall particle length ----
    ax4 = nexttile(tl);
    histogram(ax4, total_length, nbins, 'FaceColor', col.aqua, 'FaceAlpha', 0.85, ...
        'EdgeColor', edgeCol, 'LineWidth', 0.5);
    style_axes(ax4, gridCol);
    xlabel(ax4, 'total length, d_{core}+d_{bleb}-overlap (nm)'); ylabel(ax4, 'count');
    title(ax4, sprintf('Overall particle length (n = %d, %.1f\\pm%.1f nm)', ...
        numel(total_length), mean(total_length), std(total_length)));

    % ---- (5) aspect ratio ----
    ax5 = nexttile(tl);
    histogram(ax5, aspect_ratio, nbins, 'FaceColor', col.yellow, 'FaceAlpha', 0.85, ...
        'EdgeColor', edgeCol, 'LineWidth', 0.5);
    style_axes(ax5, gridCol);
    xlabel(ax5, 'aspect ratio, length / max(d_{core},d_{bleb})'); ylabel(ax5, 'count');
    title(ax5, sprintf('Aspect ratio (n = %d, %.2f\\pm%.2f)', ...
        numel(aspect_ratio), mean(aspect_ratio), std(aspect_ratio)));

    % ---- (6) bleb tilt angle, with the isotropic-sampling reference density overlaid ----
    ax6 = nexttile(tl);
    histogram(ax6, theta_tilt_deg, nbins, 'Normalization', 'pdf', 'FaceColor', col.blue, ...
        'FaceAlpha', 0.85, 'EdgeColor', edgeCol, 'LineWidth', 0.5);
    hold(ax6, 'on');
    thetaRef = linspace(0, 180, 200);
    pdfRef = (sind(thetaRef) / 2) * (pi / 180);   % density in deg^-1, for theta = acosd(2u-1)
    plot(ax6, thetaRef, pdfRef, '-', 'Color', [0.34 0.34 0.32], 'LineWidth', 1.5);
    hold(ax6, 'off');
    style_axes(ax6, gridCol);
    xlim(ax6, [0 180]);
    xlabel(ax6, 'tilt angle \theta (deg)'); ylabel(ax6, 'probability density');
    title(ax6, sprintf('Bleb tilt angle (n = %d, %.1f\\pm%.1f%s; line = isotropic ref.)', ...
        numel(theta_tilt_deg), mean(theta_tilt_deg), std(theta_tilt_deg), char(176)));

    % ---- (7) bleb azimuth angle, with the uniform-sampling reference density overlaid ----
    ax7 = nexttile(tl);
    histogram(ax7, phi_azimuth_deg, nbins, 'Normalization', 'pdf', 'FaceColor', col.orange, ...
        'FaceAlpha', 0.85, 'EdgeColor', edgeCol, 'LineWidth', 0.5);
    hold(ax7, 'on');
    yline(ax7, 1 / 360, '-', 'Color', [0.34 0.34 0.32], 'LineWidth', 1.5);
    hold(ax7, 'off');
    style_axes(ax7, gridCol);
    xlim(ax7, [0 360]);
    xlabel(ax7, 'azimuth angle \phi (deg)'); ylabel(ax7, 'probability density');
    title(ax7, sprintf('Bleb azimuth angle (n = %d, circ. mean = %.0f%s, circ. std = %.0f%s)', ...
        numel(phi_azimuth_deg), circMeanAzimuth, char(176), circStdAzimuth, char(176)));

    % ---- (8) text summary ----
    ax8 = nexttile(tl);
    axis(ax8, 'off');
    lines = {sprintf('\\bfSpherical LNPs\\rm  (n = %d)', numel(sphereLabels))};
    for g = 1:numel(groupNames)
        d_g = [sphereParams(strcmp(sphereLabels, groupNames{g})).d];
        lines{end+1} = sprintf('  %s: n=%d, %.1f\\pm%.1f nm [%.0f-%.0f]', ...
            groupNames{g}, numel(d_g), mean(d_g), std(d_g), min(d_g), max(d_g)); %#ok<AGROW>
    end
    lines = [lines, {'', sprintf('\\bfBlebbed LNPs\\rm  (n = %d)', numel(d_core)), ...
        sprintf('  core:   %.1f\\pm%.1f nm', mean(d_core), std(d_core)), ...
        sprintf('  bleb:   %.1f\\pm%.1f nm', mean(d_bleb), std(d_bleb)), ...
        sprintf('  overlap: %.1f\\pm%.1f nm', mean(neck_overlap), std(neck_overlap)), ...
        sprintf('  length: %.1f\\pm%.1f nm', mean(total_length), std(total_length)), ...
        sprintf('  aspect ratio: %.2f\\pm%.2f', mean(aspect_ratio), std(aspect_ratio)), ...
        sprintf('  tilt \\theta:   %.1f\\pm%.1f%s  (linear, isotropic ref.)', ...
            mean(theta_tilt_deg), std(theta_tilt_deg), char(176)), ...
        sprintf('  azimuth \\phi: %.0f\\pm%.0f%s  (circular mean\\pmstd)', ...
            circMeanAzimuth, circStdAzimuth, char(176))}];
    text(ax8, 0.02, 0.95, lines, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontName', 'Helvetica', 'FontSize', 10, 'Interpreter', 'tex');

    if ~isempty(p.Results.SaveAs)
        exportgraphics(fig, p.Results.SaveAs, 'Resolution', 200);
    end

    % ---- output struct ----
    stats = struct();
    for g = 1:numel(groupNames)
        d_g = [sphereParams(strcmp(sphereLabels, groupNames{g})).d];
        stats.sphere.(groupNames{g}) = struct('diameter_nm', d_g, ...
            'mean', mean(d_g), 'std', std(d_g), 'n', numel(d_g));
    end
    stats.bleb = struct('d_core_nm', d_core, 'd_bleb_nm', d_bleb, ...
        'neck_overlap_nm', neck_overlap, 'total_length_nm', total_length, ...
        'aspect_ratio', aspect_ratio, ...
        'theta_tilt_deg', theta_tilt_deg, 'theta_tilt_mean', mean(theta_tilt_deg), ...
        'theta_tilt_std', std(theta_tilt_deg), ...
        'phi_azimuth_deg', phi_azimuth_deg, ...
        'phi_azimuth_circ_mean', circMeanAzimuth, 'phi_azimuth_circ_std', circStdAzimuth, ...
        'n', numel(d_core));
    stats.fig = fig;
end

function params = collect_mesh_params(dirs, pattern)
% Read mesh.params from every file matching `pattern` in one or more folders.
    dirs = cellstr_of(dirs);
    params = struct([]);
    for i = 1:numel(dirs)
        files = dir(fullfile(dirs{i}, pattern));
        for f = 1:numel(files)
            fn = fullfile(dirs{i}, files(f).name);
            S = load(fn, 'mesh');
            if ~isfield(S, 'mesh') || ~isfield(S.mesh, 'params'), continue; end
            if isempty(params)
                params = S.mesh.params;
            else
                params(end+1) = S.mesh.params; %#ok<AGROW>
            end
        end
    end
end

function c = cellstr_of(x)
    if ischar(x) || isstring(x)
        c = {char(x)};
    else
        c = cellstr(x);
    end
end

function style_axes(ax, gridCol)
    box(ax, 'off');
    ax.XColor = [0.34 0.34 0.32];
    ax.YColor = [0.34 0.34 0.32];
    ax.GridColor = gridCol;
    ax.GridAlpha = 1;
    grid(ax, 'on');
    ax.Layer = 'top';
end

function rgb = hex2rgb(h)
    h = char(h);
    if h(1) == '#', h = h(2:end); end
    rgb = [hex2dec(h(1:2)), hex2dec(h(3:4)), hex2dec(h(5:6))] / 255;
end
