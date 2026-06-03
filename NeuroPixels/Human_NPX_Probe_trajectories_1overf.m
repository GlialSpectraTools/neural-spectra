%% probe_trajectories.m
% Plot Neuropixels probe entry points on FreeSurfer MNI pial surface
% Colored by glioma subtype: Oligo (Black), Astro (Pink), GBM (Teal)
%
% Uses FreeSurfer cvs_avg35_inMNI152 template (same as 1_ERP_results_BrainPlotting.m)

close all; clc; clear;

%% PATHS
addpath(genpath('/Users/abrahamdada/Downloads/fieldtrip-master'))
addpath('/Applications/freesurfer/8.1.0/matlab')
freesurfer_path = '/Applications/freesurfer/8.1.0/subjects/';
FS_Subject      = 'cvs_avg35_inMNI152';
surf_type       = 'pial';

output_dir = fullfile(pwd, 'plots');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% COORDINATES (from np_meta_sheet "Coords" column)
% Only recordings included in the laminar 1/f analysis

pts = {
    'NP32',  'B2',  [ 63.30,  27.00,  14.80],  'Oligo',  'R';
    'NP32',  'B3',  [ 62.40,  26.30,  17.60],  'Oligo',  'R';
    'NP37',  'B1',  [-67.48, -21.09, -25.84],  'Astro',  'L';
    'NP38',  'B5',  [-61.40,  24.20,  -0.60],  'Oligo',  'L';
    'NP38',  'B6',  [-60.50,  21.20, -10.80],  'Oligo',  'L';
    'NP46',  'B1',  [-51.70,  10.40,  33.70],  'GBM',    'L';
    'NP46',  'B2',  [-49.97,  30.54,  24.64],  'GBM',    'L';
    'NP54',  'B1',  [-61.00,  20.20,   2.40],  'GBM',    'L';
    'NP64',  'B1',  [-64.41, -28.50,   8.00],  'GBM',    'L';
    'NP65',  'B2',  [-65.40, -24.78,  -0.20],  'GBM',    'L';
    'NP78',  'B1',  [-54.80,  34.70, -16.00],  'Astro',  'L';
    'NP93',  'B1',  [-55.67,  27.68, -28.96],  'Astro',  'L';
    'NP94',  'B1',  [-64.60, -11.69,  19.27],  'GBM',    'L';
    'NP94',  'B3',  [-63.62, -26.02,  17.19],  'GBM',    'L';
    'NP95',  'B1',  [-23.98,  23.72,  51.75],  'Astro',  'L';
    'NP123', 'B1',  [-41.14,  33.49,  34.74],  'Astro',  'L';
    'NP130', 'B2',  [-60.27, -35.03,  18.97],  'GBM',    'L';
};

%% SUBTYPE COLORS (match bar plots)
clr_oligo = [0.00, 0.00, 0.00];   % Black
clr_astro = [0.91, 0.12, 0.55];   % Pink/Magenta (#E91E8C)
clr_gbm   = [0.00, 0.66, 0.62];   % Teal (#00A99D)

color_map = containers.Map();
color_map('Oligo') = clr_oligo;
color_map('Astro') = clr_astro;
color_map('GBM')   = clr_gbm;

%% LOAD FREESURFER PIAL SURFACES
SubjPath = fullfile(freesurfer_path, FS_Subject, 'surf');
% Try FieldTrip's wrapper first, fall back to FreeSurfer's read_surf
if exist('freesurfer_read_surf', 'file')
    read_fn = @freesurfer_read_surf;
elseif exist('read_surf', 'file')
    read_fn = @read_surf;
else
    error('Cannot find freesurfer_read_surf or read_surf. Check FieldTrip and FreeSurfer paths.');
end
[vtx_lh, faces_lh] = read_fn(fullfile(SubjPath, ['lh.' surf_type]));
[vtx_rh, faces_rh] = read_fn(fullfile(SubjPath, ['rh.' surf_type]));
fprintf('LH: %d vertices | RH: %d vertices\n', size(vtx_lh,1), size(vtx_rh,1));

%% ========================================================================
%  HELPER: plot brain + probes
%  ========================================================================
function plot_brain_with_probes(vtx_lh, faces_lh, pts, color_map, add_labels)
    hold on;

    bp = patch('Vertices', vtx_lh, 'Faces', faces_lh(:,[1 3 2]), ...
        'FaceColor', [0.85 0.85 0.85], 'FaceAlpha', 1, ...
        'EdgeColor', 'none');
    bp.FaceLighting    = 'gouraud';
    bp.AmbientStrength = 0.4;
    bp.DiffuseStrength = 0.7;
    material(bp, 'dull');

    for i = 1:size(pts, 1)
        subj    = pts{i, 1};
        block   = pts{i, 2};
        coord   = pts{i, 3};
        subtype = pts{i, 4};
        hemi    = pts{i, 5};

        % Mirror RH onto LH (coords are already in surface RAS)
        coord(1) = -abs(coord(1));
        pt = coord;

        % Snap to nearest vertex, then push outward along surface normal
        dists = sum((vtx_lh - pt).^2, 2);
        [~, idx] = min(dists);
        nearest = vtx_lh(idx, :);

        % Compute outward normal from neighboring faces
        face_mask = any(faces_lh == idx, 2);
        adj_faces = faces_lh(face_mask, :);
        nrm = [0 0 0];
        for f = 1:size(adj_faces, 1)
            v1 = vtx_lh(adj_faces(f,1), :);
            v2 = vtx_lh(adj_faces(f,2), :);
            v3 = vtx_lh(adj_faces(f,3), :);
            nrm = nrm + cross(v2 - v1, v3 - v1);
        end
        nrm = nrm / (norm(nrm) + eps);

        pc = nearest + 2.0 * nrm;  % 2mm outward offset

        clr = color_map(subtype);
        scatter3(pc(1), pc(2), pc(3), 130, clr, 'o', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);

        if add_labels
            lbl = sprintf('%s %s', subj, block);
            if strcmp(hemi, 'R'), lbl = [lbl ' (RH)']; end
            text(pc(1) - 2, pc(2) + 3, pc(3), lbl, ...
                'FontSize', 6.5, 'FontWeight', 'bold', 'Color', clr, ...
                'HorizontalAlignment', 'right');
        end
    end

    view(270, 0);
    camlight('headlight', 'infinite');
    lighting gouraud;
    axis equal; axis tight; axis off;
end

%% ========================================================================
%  HELPER: add manual legend
%  ========================================================================
function add_legend(clrs, names, font_size, marker_size)
    ax_leg = axes('Position', [0.65 0.05 0.30 0.14], 'Color', 'none');
    hold(ax_leg, 'on');
    for s = 1:numel(names)
        plot(ax_leg, 0.02, 1 - (s-1)*0.35, 'o', 'MarkerSize', marker_size, ...
            'MarkerFaceColor', clrs{s}, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
        text(0.08, 1 - (s-1)*0.35, names{s}, ...
            'FontSize', font_size, 'VerticalAlignment', 'middle', 'Parent', ax_leg);
    end
    set(ax_leg, 'XLim', [-0.02 0.55], 'YLim', [0 1.2]);
    axis(ax_leg, 'off');
end

%% ========================================================================
%  FIGURE 1: LABELED VERSION
%  ========================================================================
names = {'Oligodendroglioma', 'Astrocytoma', 'GBM'};
clrs  = {clr_oligo, clr_astro, clr_gbm};

fig = figure('Position', [50 50 800 700], 'Color', 'w');
plot_brain_with_probes(vtx_lh, faces_lh, pts, color_map, true);
add_legend(clrs, names, 10, 10);

print(fig, fullfile(output_dir, 'probe_trajectories_labeled.png'), '-dpng', '-r300');
print(fig, fullfile(output_dir, 'probe_trajectories_labeled.svg'), '-dsvg');
print(fig, fullfile(output_dir, 'probe_trajectories_labeled.pdf'), '-dpdf', '-bestfit');
fprintf('Saved labeled version to: %s\n', output_dir);

%% ========================================================================
%  FIGURE 2: CLEAN VERSION (publication style)
%  ========================================================================
fig2 = figure('Position', [50 50 800 700], 'Color', 'w');
plot_brain_with_probes(vtx_lh, faces_lh, pts, color_map, false);
add_legend(clrs, names, 11, 12);

print(fig2, fullfile(output_dir, 'probe_trajectories_clean.png'), '-dpng', '-r300');
print(fig2, fullfile(output_dir, 'probe_trajectories_clean.svg'), '-dsvg');
print(fig2, fullfile(output_dir, 'probe_trajectories_clean.pdf'), '-dpdf', '-bestfit');
fprintf('Saved clean version to: %s\n', output_dir);

%% ========================================================================
%  FIGURE CAPTION
%  ========================================================================
caption = [ ...
    'Figure X. Neuropixels recording locations across glioma-infiltrated cortex. ' ...
    'Cortical entry points of Neuropixels 1.0 probes plotted on a FreeSurfer MNI152 ' ...
    'template pial surface (left lateral view; right-hemisphere recordings mirrored). ' ...
    'Dots are colored by glioma subtype: oligodendroglioma (black), astrocytoma (pink), ' ...
    'and glioblastoma (teal). Recordings were obtained from lateral frontal, temporal, ' ...
    'and parietal cortices in n = 13 patients (17 unique trajectories, 21 probe insertions) ' ...
    'with WHO grade 2-4 diffuse glioma, including oligodendroglioma (O; n = 2 patients, ' ...
    '4 trajectories, 4 insertions), astrocytoma (A; n = 5 patients, 5 trajectories, ' ...
    '6 insertions), and glioblastoma (GBM; n = 6 patients, 8 trajectories, 11 insertions). ' ...
    'Each dot represents one unique cortical entry point; fewer dots than total insertions ' ...
    'reflect multi-probe or repeated penetrations at the same site.'];

fid = fopen(fullfile(output_dir, 'probe_trajectories_caption.txt'), 'w');
fprintf(fid, '%s\n', caption);
fclose(fid);
fprintf('Saved caption to: %s\n', fullfile(output_dir, 'probe_trajectories_caption.txt'));
