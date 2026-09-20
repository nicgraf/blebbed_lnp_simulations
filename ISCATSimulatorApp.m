classdef ISCATSimulatorApp < matlab.apps.AppBase
    % ISCATSIMULATORAPP  Point-and-click iSCAT simulator app.
    %
    % Run it just like any MATLAB app: open this file and click Run (or
    % type ISCATSimulatorApp at the command line). A window opens where
    % you can:
    %   1. Point it at your "lnp_meshes" folder and load + solve the
    %      particle meshes (same physics/specs as simulation.mlx).
    %   2. Pick which particle to look at (population + particle number).
    %   3. Click Simulate -- it auto-focuses and renders both the
    %      in-plane (best-focus) image AND the full z-stack (a 2D image
    %      at every defocus step, not just a 1D profile), which you can
    %      scrub through with the slider.
    %   4. Download the in-plane image as a PNG, and the z-stack as a
    %      multi-page TIFF (viewable in ImageJ/Fiji) plus a .mat file
    %      with the raw numeric data, via the Save buttons.

    properties (Access = public)
        UIFigure              matlab.ui.Figure
        MainGrid              matlab.ui.container.GridLayout
        LeftPanel              matlab.ui.container.Panel
        LeftGrid               matlab.ui.container.GridLayout
        TitleLabel               matlab.ui.control.Label
        MeshFolderLabel          matlab.ui.control.Label
        MeshFolderRow            matlab.ui.container.GridLayout
        MeshFolderEditField        matlab.ui.control.EditField
        BrowseButton                matlab.ui.control.Button
        LoadMeshesButton         matlab.ui.control.Button
        PopulationLabel          matlab.ui.control.Label
        PopulationDropdown       matlab.ui.control.DropDown
        ParticleLabel            matlab.ui.control.Label
        ParticleSpinner          matlab.ui.control.Spinner
        SimulateButton           matlab.ui.control.Button
        ExportLabel              matlab.ui.control.Label
        ExportImageButton        matlab.ui.control.Button
        ExportStackButton        matlab.ui.control.Button
        LogLabel                 matlab.ui.control.Label
        StatusTextArea           matlab.ui.control.TextArea
        RightGrid              matlab.ui.container.GridLayout
        InPlaneAxes              matlab.ui.control.UIAxes
        ZStackAxes               matlab.ui.control.UIAxes
        SliderGrid              matlab.ui.container.GridLayout
        ZSlider                    matlab.ui.control.Slider
        ZValueLabel                matlab.ui.control.Label
    end

    properties (Access = private)
        full_lnp = {}
        empty_lnp = {}
        blebbed_lnp = {}
        specs

        stack = []          % npix x npix x n_z contrast volume, last simulated particle
        stack_clims = [-1 1]
        x_nm = []
        z_nm = []
        best_z = 0
        best_z_index = 1
        current_tag = ''
    end

    methods (Access = private)

        %% ---- small logging helper: writes to the on-screen log AND the command window ----
        function log(app, msg)
            fprintf('%s\n', msg);
            existing = app.StatusTextArea.Value;
            if ischar(existing)
                existing = {existing};
            end
            app.StatusTextArea.Value = [existing; {msg}];
            drawnow limitrate;
        end

        function StartupFcn(app)
            addpath(fileparts(mfilename('fullpath')));   % so bleb_volume.m is found
            app.specs = default_specs();

            default_dir = fullfile(fileparts(mfilename('fullpath')), 'lnp_meshes');
            if exist(default_dir, 'dir')
                app.MeshFolderEditField.Value = default_dir;
            end

            if exist('Material', 'file') ~= 2 && exist('Material', 'class') ~= 8
                uialert(app.UIFigure, ...
                    ['The "nanobem" toolbox does not appear to be on your MATLAB path. ' ...
                     'Ask whoever set up your MATLAB for its location, then add it via ' ...
                     'Home tab > Set Path > Add with Subfolders, then run this app again.'], ...
                    'nanobem not found');
            end

            app.log('Ready. Choose a mesh folder and click "Load & Solve Meshes".');
        end

        function BrowseButtonPushed(app, ~)
            folder = uigetdir(app.MeshFolderEditField.Value, 'Select the LNP mesh folder');
            if ~isequal(folder, 0)
                app.MeshFolderEditField.Value = folder;
            end
        end

        function LoadMeshesButtonPushed(app, ~)
            mesh_dir = app.MeshFolderEditField.Value;
            if isempty(mesh_dir) || ~exist(mesh_dir, 'dir')
                uialert(app.UIFigure, 'Please choose a valid mesh folder first.', 'No mesh folder');
                return;
            end

            d = uiprogressdlg(app.UIFigure, 'Title', 'Loading meshes', ...
                'Message', 'Loading mesh files...', 'Indeterminate', 'on', 'Cancelable', false);
            cleanupD = onCleanup(@() close(d));

            try
                app.log(sprintf('Loading meshes from: %s', mesh_dir));
                app.full_lnp = normalize_population( ...
                    importdata(pick_mesh_file(mesh_dir, 'full_lnp')), {'tau', 'diameter'}, 'full_lnp');
                app.empty_lnp = normalize_population( ...
                    importdata(pick_mesh_file(mesh_dir, 'empty_lnp')), {'tau', 'diameter'}, 'empty_lnp');
                app.blebbed_lnp = normalize_population( ...
                    importdata(pick_mesh_file(mesh_dir, 'blebbed_lnp')), ...
                    {'tau_lipid_out', 'tau_aqueous_out', 'tau_interface', 'd_core', 'd_bleb', 'neck_overlap', ...
                     'phi_azimuth_deg', 'theta_tilt_deg'}, 'blebbed_lnp');
                app.log(sprintf('  full_lnp: %d, empty_lnp: %d, blebbed_lnp: %d particles', ...
                    numel(app.full_lnp), numel(app.empty_lnp), numel(app.blebbed_lnp)));

                d.Indeterminate = 'off';
                n_total = numel(app.full_lnp) + numel(app.empty_lnp) + numel(app.blebbed_lnp);
                n_done = 0;

                for i = 1:numel(app.full_lnp)
                    tau  = app.full_lnp{i}.tau;
                    qinc = app.specs.einc(tau, 'layer', app.specs.layer);
                    bem  = stratified.bemsolver(tau, app.specs.layer, 'order', []);
                    app.full_lnp{i}.sol = bem \ qinc;
                    n_done = n_done + 1;
                    d.Value = n_done / n_total;
                    d.Message = sprintf('Solving full lnp (%d/%d)...', i, numel(app.full_lnp));
                end

                for i = 1:numel(app.empty_lnp)
                    tau  = app.empty_lnp{i}.tau;
                    qinc = app.specs.einc(tau, 'layer', app.specs.layer);
                    bem  = stratified.bemsolver(tau, app.specs.layer, 'order', []);
                    app.empty_lnp{i}.sol = bem \ qinc;
                    n_done = n_done + 1;
                    d.Value = n_done / n_total;
                    d.Message = sprintf('Solving empty lnp (%d/%d)...', i, numel(app.empty_lnp));
                end

                for i = 1:numel(app.blebbed_lnp)
                    tau = merge_bleb_patches(app.blebbed_lnp{i}.tau_lipid_out, ...
                                             app.blebbed_lnp{i}.tau_aqueous_out, ...
                                             app.blebbed_lnp{i}.tau_interface);
                    qinc = app.specs.einc(tau, 'layer', app.specs.layer);
                    bem  = stratified.bemsolver(tau, app.specs.layer, 'order', 5);
                    app.blebbed_lnp{i}.tau = tau;
                    app.blebbed_lnp{i}.sol = bem \ qinc;
                    n_done = n_done + 1;
                    d.Value = n_done / n_total;
                    d.Message = sprintf('Solving blebbed lnp (%d/%d)...', i, numel(app.blebbed_lnp));
                end
            catch ME
                app.log(['ERROR: ' ME.message]);
                uialert(app.UIFigure, ME.message, 'Failed to load/solve meshes');
                return;
            end

            app.log('BEM solve complete. Ready to simulate.');
            app.PopulationDropdown.Enable = 'on';
            app.ParticleSpinner.Enable = 'on';
            app.SimulateButton.Enable = 'on';
            app.PopulationDropdownValueChanged();
        end

        function PopulationDropdownValueChanged(app, ~)
            n = numel(app.currentPopData());
            if n < 1
                n = 1;
            end
            app.ParticleSpinner.Limits = [1, n];
            app.ParticleSpinner.Value = min(app.ParticleSpinner.Value, n);
        end

        function pop_data = currentPopData(app)
            switch app.PopulationDropdown.Value
                case 'Full LNP'
                    pop_data = app.full_lnp;
                case 'Empty LNP'
                    pop_data = app.empty_lnp;
                otherwise
                    pop_data = app.blebbed_lnp;
            end
        end

        function SimulateButtonPushed(app, ~)
            pop_data = app.currentPopData();
            if isempty(pop_data)
                uialert(app.UIFigure, 'Load meshes first.', 'Nothing to simulate');
                return;
            end
            idx = round(app.ParticleSpinner.Value);
            particle = pop_data{idx};
            pop_name = app.PopulationDropdown.Value;

            d = uiprogressdlg(app.UIFigure, 'Title', 'Simulating', ...
                'Message', 'Computing far fields...', 'Value', 0, 'Cancelable', false);
            cleanupD = onCleanup(@() close(d));

            try
                [stack_out, x, z_scan, bz] = render_iscat_zstack(particle.sol, app.specs, d);
            catch ME
                app.log(['ERROR: ' ME.message]);
                uialert(app.UIFigure, ME.message, 'Simulation failed');
                return;
            end

            app.stack = stack_out;
            app.x_nm = x;
            app.z_nm = z_scan;
            app.best_z = bz;
            [~, app.best_z_index] = min(abs(z_scan - bz));
            app.current_tag = sprintf('%s_%d', strrep(pop_name, ' ', ''), idx);

            cmax = max(abs(app.stack(:)));
            if cmax == 0
                cmax = 1;
            end
            app.stack_clims = [-cmax, cmax];

            if isfield(particle, 'diameter')
                size_label = sprintf('diameter=%.0f nm', particle.diameter);
            else
                d_eq = (6*bleb_volume(particle.d_core, particle.d_bleb, particle.neck_overlap)/pi)^(1/3);
                size_label = sprintf('d_{eq}=%.0f nm, az=%.0f%s, tilt=%.0f%s', ...
                    d_eq, particle.phi_azimuth_deg, char(176), particle.theta_tilt_deg, char(176));
            end

            imagesc(app.InPlaneAxes, app.x_nm, app.x_nm, app.stack(:, :, app.best_z_index), app.stack_clims);
            axis(app.InPlaneAxes, 'image');
            colormap(app.InPlaneAxes, 'gray');
            colorbar(app.InPlaneAxes);
            xlabel(app.InPlaneAxes, 'x (nm)'); ylabel(app.InPlaneAxes, 'y (nm)');
            title(app.InPlaneAxes, sprintf('%s #%d, %s, z=%.0f nm', pop_name, idx, size_label, app.best_z));

            app.ZSlider.Limits = [1, numel(app.z_nm)];
            app.ZSlider.Value = app.best_z_index;
            app.ZSlider.MajorTicks = round(linspace(1, numel(app.z_nm), 5));
            app.ZSliderValueChanged([]);

            app.ExportImageButton.Enable = 'on';
            app.ExportStackButton.Enable = 'on';

            app.log(sprintf('Simulated %s #%d -- best focus z = %.0f nm. Full z-stack: %d planes, %.2f um deep total.', ...
                pop_name, idx, app.best_z, numel(app.z_nm), (max(app.z_nm)-min(app.z_nm))/1000));
        end

        function ZSliderValueChanged(app, event)
            if isempty(app.stack)
                return;
            end
            if nargin >= 2 && ~isempty(event) && isprop(event, 'Value')
                iz = round(event.Value);   % live value while dragging (ValueChangingFcn)
            else
                iz = round(app.ZSlider.Value);
            end
            iz = max(1, min(numel(app.z_nm), iz));
            imagesc(app.ZStackAxes, app.x_nm, app.x_nm, app.stack(:, :, iz), app.stack_clims);
            axis(app.ZStackAxes, 'image');
            colormap(app.ZStackAxes, 'gray');
            colorbar(app.ZStackAxes);
            xlabel(app.ZStackAxes, 'x (nm)'); ylabel(app.ZStackAxes, 'y (nm)');
            tag = '';
            if iz == app.best_z_index
                tag = '  (best focus)';
            end
            title(app.ZStackAxes, sprintf('z = %.0f nm%s', app.z_nm(iz), tag));
            app.ZValueLabel.Text = sprintf('z = %.0f nm  (plane %d of %d)', app.z_nm(iz), iz, numel(app.z_nm));
        end

        function ExportImageButtonPushed(app, ~)
            if isempty(app.stack)
                return;
            end
            [file, path] = uiputfile({'*.png', 'PNG image'}, 'Save in-plane image', ...
                sprintf('iSCAT_%s_inplane.png', app.current_tag));
            if isequal(file, 0)
                return;
            end
            exportgraphics(app.InPlaneAxes, fullfile(path, file), 'Resolution', 200);
            app.log(['Saved in-plane image: ' fullfile(path, file)]);
        end

        function ExportStackButtonPushed(app, ~)
            if isempty(app.stack)
                return;
            end
            [file, path] = uiputfile({'*.tif', 'Multi-page TIFF (z-stack, for ImageJ/Fiji)'}, ...
                'Save z-stack', sprintf('iSCAT_%s_zstack.tif', app.current_tag));
            if isequal(file, 0)
                return;
            end
            [~, base_name] = fileparts(file);
            tif_path = fullfile(path, [base_name '.tif']);
            mat_path = fullfile(path, [base_name '.mat']);

            cmax = max(abs(app.stack(:)));
            if cmax == 0
                cmax = 1;
            end
            scaled = uint16((app.stack + cmax) / (2 * cmax) * 65535);
            for iz = 1:size(scaled, 3)
                if iz == 1
                    imwrite(scaled(:, :, iz), tif_path);
                else
                    imwrite(scaled(:, :, iz), tif_path, 'WriteMode', 'append');
                end
            end

            contrast_stack_percent = app.stack;   %#ok<NASGU> raw iSCAT contrast, percent, npix x npix x n_z
            x_nm = app.x_nm;                       %#ok<NASGU>
            z_nm = app.z_nm;                       %#ok<NASGU>
            best_z_nm = app.best_z;                %#ok<NASGU>
            tag = app.current_tag;                 %#ok<NASGU>
            save(mat_path, 'contrast_stack_percent', 'x_nm', 'z_nm', 'best_z_nm', 'tag');

            app.log(sprintf('Saved z-stack: %s (viewable TIFF) and %s (raw data)', tif_path, mat_path));
        end

        function UIFigureCloseRequest(app, ~)
            delete(app);
        end

        %% ---- UI layout ----
        function createComponents(app)
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1100 680];
            app.UIFigure.Name = 'iSCAT Simulator';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            app.MainGrid = uigridlayout(app.UIFigure, [1, 2]);
            app.MainGrid.ColumnWidth = {320, '1x'};
            app.MainGrid.Padding = [10 10 10 10];

            %% ---- left control panel ----
            app.LeftPanel = uipanel(app.MainGrid);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            app.LeftGrid = uigridlayout(app.LeftPanel);
            app.LeftGrid.ColumnWidth = {'1x'};
            app.LeftGrid.RowHeight = {30, 20, 30, 32, 12, 20, 30, 20, 30, 32, 12, 20, 30, 30, 12, 20, '1x'};
            app.LeftGrid.RowSpacing = 4;

            app.TitleLabel = uilabel(app.LeftGrid);
            app.TitleLabel.Text = 'iSCAT Simulator';
            app.TitleLabel.FontSize = 18;
            app.TitleLabel.FontWeight = 'bold';
            app.TitleLabel.Layout.Row = 1; app.TitleLabel.Layout.Column = 1;

            app.MeshFolderLabel = uilabel(app.LeftGrid);
            app.MeshFolderLabel.Text = 'Mesh folder:';
            app.MeshFolderLabel.Layout.Row = 2; app.MeshFolderLabel.Layout.Column = 1;

            app.MeshFolderRow = uigridlayout(app.LeftGrid, [1, 2]);
            app.MeshFolderRow.ColumnWidth = {'1x', 70};
            app.MeshFolderRow.Padding = [0 0 0 0];
            app.MeshFolderRow.Layout.Row = 3; app.MeshFolderRow.Layout.Column = 1;

            app.MeshFolderEditField = uieditfield(app.MeshFolderRow, 'text');
            app.MeshFolderEditField.Layout.Row = 1; app.MeshFolderEditField.Layout.Column = 1;

            app.BrowseButton = uibutton(app.MeshFolderRow, 'push');
            app.BrowseButton.Text = 'Browse...';
            app.BrowseButton.Layout.Row = 1; app.BrowseButton.Layout.Column = 2;
            app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);

            app.LoadMeshesButton = uibutton(app.LeftGrid, 'push');
            app.LoadMeshesButton.Text = 'Load & Solve Meshes';
            app.LoadMeshesButton.Layout.Row = 4; app.LoadMeshesButton.Layout.Column = 1;
            app.LoadMeshesButton.ButtonPushedFcn = createCallbackFcn(app, @LoadMeshesButtonPushed, true);

            app.PopulationLabel = uilabel(app.LeftGrid);
            app.PopulationLabel.Text = 'Population:';
            app.PopulationLabel.Layout.Row = 6; app.PopulationLabel.Layout.Column = 1;

            app.PopulationDropdown = uidropdown(app.LeftGrid);
            app.PopulationDropdown.Items = {'Full LNP', 'Empty LNP', 'Blebbed LNP'};
            app.PopulationDropdown.Enable = 'off';
            app.PopulationDropdown.Layout.Row = 7; app.PopulationDropdown.Layout.Column = 1;
            app.PopulationDropdown.ValueChangedFcn = createCallbackFcn(app, @PopulationDropdownValueChanged, true);

            app.ParticleLabel = uilabel(app.LeftGrid);
            app.ParticleLabel.Text = 'Particle #:';
            app.ParticleLabel.Layout.Row = 8; app.ParticleLabel.Layout.Column = 1;

            app.ParticleSpinner = uispinner(app.LeftGrid);
            app.ParticleSpinner.Limits = [1 1];
            app.ParticleSpinner.RoundFractionalValues = 'on';
            app.ParticleSpinner.Step = 1;
            app.ParticleSpinner.Enable = 'off';
            app.ParticleSpinner.Layout.Row = 9; app.ParticleSpinner.Layout.Column = 1;

            app.SimulateButton = uibutton(app.LeftGrid, 'push');
            app.SimulateButton.Text = 'Simulate';
            app.SimulateButton.FontWeight = 'bold';
            app.SimulateButton.Enable = 'off';
            app.SimulateButton.Layout.Row = 10; app.SimulateButton.Layout.Column = 1;
            app.SimulateButton.ButtonPushedFcn = createCallbackFcn(app, @SimulateButtonPushed, true);

            app.ExportLabel = uilabel(app.LeftGrid);
            app.ExportLabel.Text = 'Download results:';
            app.ExportLabel.Layout.Row = 12; app.ExportLabel.Layout.Column = 1;

            app.ExportImageButton = uibutton(app.LeftGrid, 'push');
            app.ExportImageButton.Text = 'Save In-Plane Image (PNG)';
            app.ExportImageButton.Enable = 'off';
            app.ExportImageButton.Layout.Row = 13; app.ExportImageButton.Layout.Column = 1;
            app.ExportImageButton.ButtonPushedFcn = createCallbackFcn(app, @ExportImageButtonPushed, true);

            app.ExportStackButton = uibutton(app.LeftGrid, 'push');
            app.ExportStackButton.Text = 'Save Z-Stack (TIFF + MAT)';
            app.ExportStackButton.Enable = 'off';
            app.ExportStackButton.Layout.Row = 14; app.ExportStackButton.Layout.Column = 1;
            app.ExportStackButton.ButtonPushedFcn = createCallbackFcn(app, @ExportStackButtonPushed, true);

            app.LogLabel = uilabel(app.LeftGrid);
            app.LogLabel.Text = 'Log:';
            app.LogLabel.Layout.Row = 16; app.LogLabel.Layout.Column = 1;

            app.StatusTextArea = uitextarea(app.LeftGrid);
            app.StatusTextArea.Editable = 'off';
            app.StatusTextArea.Layout.Row = 17; app.StatusTextArea.Layout.Column = 1;

            %% ---- right display area ----
            app.RightGrid = uigridlayout(app.MainGrid, [2, 2]);
            app.RightGrid.Layout.Row = 1; app.RightGrid.Layout.Column = 2;
            app.RightGrid.RowHeight = {'1x', 60};
            app.RightGrid.ColumnWidth = {'1x', '1x'};

            app.InPlaneAxes = uiaxes(app.RightGrid);
            app.InPlaneAxes.Layout.Row = 1; app.InPlaneAxes.Layout.Column = 1;
            title(app.InPlaneAxes, 'In-plane image (best focus)');

            app.ZStackAxes = uiaxes(app.RightGrid);
            app.ZStackAxes.Layout.Row = 1; app.ZStackAxes.Layout.Column = 2;
            title(app.ZStackAxes, 'Z-stack viewer (drag slider below)');

            app.SliderGrid = uigridlayout(app.RightGrid, [1, 2]);
            app.SliderGrid.Layout.Row = 2; app.SliderGrid.Layout.Column = [1 2];
            app.SliderGrid.ColumnWidth = {'1x', 220};
            app.SliderGrid.Padding = [10 10 10 0];

            app.ZSlider = uislider(app.SliderGrid);
            app.ZSlider.Limits = [1 2];
            app.ZSlider.Layout.Row = 1; app.ZSlider.Layout.Column = 1;
            app.ZSlider.ValueChangedFcn = createCallbackFcn(app, @ZSliderValueChanged, true);
            app.ZSlider.ValueChangingFcn = createCallbackFcn(app, @ZSliderValueChanged, true);

            app.ZValueLabel = uilabel(app.SliderGrid);
            app.ZValueLabel.Text = 'z = -- nm';
            app.ZValueLabel.Layout.Row = 1; app.ZValueLabel.Layout.Column = 2;

            app.UIFigure.Visible = 'on';
        end
    end

    methods (Access = public)
        function app = ISCATSimulatorApp()
            createComponents(app)
            registerApp(app, app.UIFigure)
            runStartupFcn(app, @StartupFcn)
            if nargout == 0
                clear app
            end
        end

        function delete(app)
            delete(app.UIFigure)
        end
    end
end


%% ======================================================================
%  Local functions (physics/specs -- identical to simulation.mlx)
%  ======================================================================

function specs = default_specs()
% All the physical/optical settings below match the lab's real microscope
% and the original simulation pipeline exactly (simulation.mlx).

    specs.lambda = 555;              % nm
    specs.k0 = 2 * pi / specs.lambda;

    specs.n_glass = 1.515;
    specs.n_water = 1.335;
    specs.mat_glass = Material(specs.n_glass^2, 1);
    specs.mat_water = Material(specs.n_water^2, 1);
    specs.layer = stratified.layerstructure([specs.mat_glass, specs.mat_water], 0);

    specs.pol = [1, 0, 0];
    specs.dir = [0, 0, 1];
    specs.einc = optics.decompose(specs.k0, specs.pol, specs.dir);

    specs.NA = 1.45;
    specs.mag_objective = 100;
    specs.mag_intermediate = 1.5;
    specs.mag_total = specs.mag_objective * specs.mag_intermediate;   % 150x

    specs.sensor_fov_um = 118;
    specs.sensor_pixels = 1608;
    specs.pixel_size_nm = specs.sensor_fov_um * 1000 / specs.sensor_pixels;   % ~73.38 nm/px

    specs.npix = 61;
    specs.icenter = ceil(specs.npix / 2);
    specs.x = specs.pixel_size_nm * ( -(specs.npix-1)/2 : (specs.npix-1)/2 );   % nm

    specs.air = Material(1, 1);
    specs.rot = optics.roty(180);
    specs.lens = optics.lensimage(specs.mat_glass, specs.air, specs.k0, specs.NA, 'rot', specs.rot);
    specs.refl = secondary(specs.einc, specs.layer, 'dir', 'down');

    % Full depth of the z-stack: +/-2 um around focus (4 um total), same
    % range used for defocus scans in simulation.mlx.
    specs.z_scan = -2000:20:2000;    % nm
end

function f = pick_mesh_file(mesh_dir, base_name)
% Prefers "<base_name>_sol.mat" (matches the lab's original file naming)
% but falls back to "<base_name>.mat" if that's what's present.
    candidates = {[base_name '_sol.mat'], [base_name '.mat']};
    for i = 1:numel(candidates)
        f = fullfile(mesh_dir, candidates{i});
        if exist(f, 'file')
            return;
        end
    end
    error('ISCATSimulatorApp:meshFileNotFound', ...
        'Could not find %s_sol.mat or %s.mat in %s', base_name, base_name, mesh_dir);
end

function tau = merge_bleb_patches(tau_lipid_out, tau_aqueous_out, tau_interface)
% MERGE_BLEB_PATCHES  Combine the three separately-built blebbed-LNP
% patches into ONE BoundaryEdge with a single, consistent set of global
% edge indices.
%
% Why this matters: BoundaryEdge numbers its edge basis functions (nu)
% starting from 1 every time it is built. Three separately-built patches
% therefore all reuse indices 1..N, so simply concatenating them makes the
% solver silently add unrelated unknowns together, and the shared rim
% edges where the patches meet are never linked. Re-running nanobem's own
% ShapeEdge numbering over ALL elements at once (exactly what
% BoundaryEdge(mat, p1, inout1, p2, inout2, ...) does internally) fixes
% both. Works on the stored meshes -- no need to regenerate them.

    tau_all = [tau_lipid_out(:).', tau_aqueous_out(:).', tau_interface(:).'];
    shape = ShapeEdge(tau_all);
    tau = BoundaryEdge(shape);
    %  keep the same element order as TAU_ALL (same as BoundaryEdge/init2)
    [~, i1] = ismember(vertcat(tau_all.pos), vertcat(tau.pos), 'rows');
    tau = tau(i1);
end

function pop_data = normalize_population(raw, required_fields, label)
% NORMALIZE_POPULATION  Coerce a loaded population variable into the
% Nx1 cell-of-structs shape the rest of the app expects.
%
% Tolerates two common variants in how these .mat files can end up
% getting saved: a struct array instead of a cell array of structs, and
% each particle accidentally double-wrapped in an extra 1x1 cell (e.g.
% from a cellfun/arrayfun call that wasn't given 'UniformOutput', false
% at the right level). Raises a clear, specific error -- naming the file,
% the particle index, and what was actually found -- if the shape still
% isn't recognizable, instead of letting a cryptic "dot indexing" error
% surface deep inside the BEM solve loop.

    if isstruct(raw)
        raw = num2cell(raw(:));
    end

    if ~iscell(raw)
        error('ISCATSimulatorApp:badPopulationShape', ...
            ['%s: expected a cell array (or struct array) of particles in the ' ...
             '.mat file, but found class ''%s'' instead.'], label, class(raw));
    end

    raw = raw(:);
    pop_data = cell(numel(raw), 1);
    for i = 1:numel(raw)
        item = raw{i};
        while iscell(item) && isscalar(item)   % peel off accidental extra wrapping
            item = item{1};
        end
        if ~isstruct(item) || ~isscalar(item)
            error('ISCATSimulatorApp:badPopulationShape', ...
                ['%s: particle %d is not a single struct (found class ''%s''). ' ...
                 'Check how this .mat file was saved.'], label, i, class(item));
        end
        missing = required_fields(~isfield(item, required_fields));
        if ~isempty(missing)
            error('ISCATSimulatorApp:missingFields', ...
                '%s: particle %d is missing field(s): %s', label, i, strjoin(missing, ', '));
        end
        pop_data{i} = item;
    end
end

function [stack, x, z_scan, best_z] = render_iscat_zstack(sol, specs, progressDlg)
% RENDER_ISCAT_ZSTACK  Full 2D iSCAT contrast image at every z in
% specs.z_scan (not just the central-row profile used in
% simulation.mlx's render_iscat_zscan) -- this is the "full z-stack".
% The best-focus plane is picked with the exact same criterion as
% simulation.mlx: the z where the raw center-pixel interference signal
% peaks, excluding any z where the reference beam has collapsed too low
% to trust.

    far = farfields(sol, specs.lens.dir);
    n_z = numel(specs.z_scan);
    npix = numel(specs.x);
    stack = zeros(npix, npix, n_z);
    ibg_center = zeros(n_z, 1);
    numerator_center = zeros(n_z, 1);

    for iz = 1:n_z
        if nargin >= 3 && ~isempty(progressDlg)
            progressDlg.Value = iz / n_z;
            progressDlg.Message = sprintf('Rendering z = %.0f nm (%d/%d)...', specs.z_scan(iz), iz, n_z);
        end

        focus = [0, 0, specs.z_scan(iz)];
        isca = efield(specs.lens, far,        specs.x, specs.x, 'focus', focus);
        iref = efield(specs.lens, specs.refl, specs.x, specs.x, 'focus', focus);
        isca = flip(isca, 1);   % consequence of the rot=180 lens axis
        iref = flip(iref, 1);

        im  = dot(isca + iref, isca + iref, 3);
        ibg = dot(iref, iref, 3);
        contrast = (im - ibg) ./ ibg * 100;

        stack(:, :, iz) = contrast;
        ibg_center(iz)       = ibg(specs.icenter, specs.icenter);
        numerator_center(iz) = im(specs.icenter, specs.icenter) - ibg(specs.icenter, specs.icenter);
    end

    x = specs.x;
    z_scan = specs.z_scan;

    ibg_floor_frac = 0.1;
    ibg_ok = ibg_center >= ibg_floor_frac * max(ibg_center);
    numerator_masked = abs(numerator_center);
    numerator_masked(~ibg_ok) = -Inf;
    if all(~ibg_ok)
        numerator_masked = abs(numerator_center);
    end
    [~, idx_best] = max(numerator_masked);
    best_z = z_scan(idx_best);
end
