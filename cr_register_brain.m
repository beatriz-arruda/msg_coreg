% cr_register_brain - Register SPM template brain meshes to a subject head surface
%
% Places the SPM canonical cortex, scalp, inner skull and outer skull in
% the subject's space by matching the template fiducials (nasion, left and
% right pre-auricular points) to fiducials selected on the subject's head.
%
% The transform is a rigid-body fit of the three fiducial pairs (rotation
% and translation, least squares), optionally with one isotropic scale
% factor. It never mirrors the template: a mirrored brain would swap the
% hemispheres. Fiducials selected in the wrong order are the usual cause of
% a brain that ends up in the neck or outside the head; the two orders of
% the pre-auricular points are therefore both tried, and if only the
% swapped order puts the brain inside the head, that order is used with a
% warning.
%
% USAGE:
%   [temp_brain, T, info] = cr_register_brain(S)
%
% INPUT:
%   S              - Structure with the following fields:
%
%   Required:
%     S.fiducials    - 3x3 subject fiducial coordinates, one point per row,
%                      in the units of S.head (normally mm):
%                        Row 1: Nasion (NAS)
%                        Row 2: Left pre-auricular point (LPA)
%                        Row 3: Right pre-auricular point (RPA)
%     S.head         - Subject head/torso mesh (.vertices, .faces); used to
%                      check that the brain lands inside the head
%
%   Optional:
%     S.scale        - true: also scale the template isotropically so its
%                      fiducial triangle matches the subject's in size
%                      (default: false, rigid only)
%     S.check_order  - try LPA/RPA in both orders (default: true)
%     S.plot         - display the registration (default: false)
%
% OUTPUT:
%   temp_brain     - registered SPM canonical meshes (.vertices, .faces):
%                      .brain   cortical surface
%                      .scalp   scalp
%                      .iskull  inner skull
%                      .oskull  outer skull
%   T              - 4x4 transform, template -> subject (homogeneous, row
%                    vectors transformed as (T * [v 1]')')
%   info           - quality of the fit:
%                      .fid_error_mm   [3x1] distance between each
%                                      subject fiducial and its template
%                                      fiducial after registration
%                      .frac_inside    fraction of cortex vertices inside
%                                      S.head
%                      .swapped        true if LPA and RPA were swapped
%                      .scale          isotropic scale applied (1 if none)
%
% NOTES:
%   - Template fiducials come from SPM's EEGtemplates/fiducials.sfp.
%   - Units: the template is in mm. If the subject is in m or cm, a power
%     of ten is applied first, estimated from the fiducial triangle sizes.
%   - Fiducial errors of a few mm are normal (template and subject heads
%     differ in shape). Errors of several cm, or frac_inside well below 1,
%     mean the fiducials were misplaced; select them again.
%
% EXAMPLE:
%   S.fiducials = [nas; lpa; rpa];   % 3x3, one point per row
%   S.head      = torso_mesh;
%   S.plot      = true;
%   [temp_brain, T, info] = cr_register_brain(S);
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_coreg
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Coregistration Toolbox.

function [temp_brain, T, info] = cr_register_brain(S)

if ~isfield(S, 'fiducials'), error('You must provide fiducial locations (NAS, LPA, RPA)'); end
if ~isfield(S, 'head'), error('You must provide a head shape to fit to!'); end
if ~isfield(S, 'scale'), S.scale = false; end
if ~isfield(S, 'check_order'), S.check_order = true; end
if ~isfield(S, 'plot'), S.plot = false; end
if ~isequal(size(S.fiducials), [3 3])
    error('S.fiducials must be 3x3 (rows NAS, LPA, RPA); got %dx%d.', size(S.fiducials, 1), size(S.fiducials, 2));
end
fid_subject = double(S.fiducials);

% SPM template meshes and fiducials (mm)
meshes = struct();
files = {'brain', 'cortex_5124'; 'scalp', 'scalp_2562'; 'iskull', 'iskull_2562'; 'oskull', 'oskull_2562'};
for k = 1:size(files, 1)
    g = gifti(fullfile(spm('dir'), 'canonical', [files{k, 2} '.surf.gii']));
    meshes.(files{k, 1}) = struct('vertices', double(g.vertices), 'faces', double(g.faces));
end
fid_template = spm_eeg_fixpnt(ft_read_headshape( ...
    fullfile(spm('dir'), 'EEGtemplates', 'fiducials.sfp')));
fid_template = fid_template.fid.pnt(1:3, :);   % NAS, LPA, RPA

% Fit with the fiducials as given, and with LPA and RPA swapped
orders = {[1 2 3]};
if S.check_order, orders{end + 1} = [1 3 2]; end
fits = cell(size(orders));
for o = 1:numel(orders)
    fits{o} = fit_fiducials(fid_subject(orders{o}, :), fid_template, S.scale, ...
        meshes.brain.vertices, S.head);
end

use = 1;
if numel(fits) > 1 && fits{1}.frac_inside < 0.5 && fits{2}.frac_inside > max(0.9, fits{1}.frac_inside)
    use = 2;
    warning('cr_register_brain:order', ['LPA and RPA appear to be swapped: with the order given only ' ...
        '%.0f%% of the brain is inside the head, with them swapped %.0f%%. Using the swapped order; ' ...
        'select NAS, LPA (subject''s left), RPA (subject''s right) to avoid this.'], ...
        100 * fits{1}.frac_inside, 100 * fits{2}.frac_inside);
end
F = fits{use};
T = F.T;

info.fid_error_mm = F.fid_error;
info.frac_inside  = F.frac_inside;
info.swapped      = use == 2;
info.scale        = F.scale;
fprintf('Brain registration: fiducial error NAS %.1f, LPA %.1f, RPA %.1f mm; %.0f%% of cortex inside the head\n', ...
    F.fid_error, 100 * F.frac_inside);
if F.frac_inside < 0.9
    warning('cr_register_brain:outside', ['Only %.0f%% of the cortex is inside the head after registration. ' ...
        'Check the fiducials (order NAS, LPA, RPA; LPA on the subject''s left).'], 100 * F.frac_inside);
end
if max(F.fid_error) > 20
    warning('cr_register_brain:fiderror', ['Fiducial error up to %.0f mm: the subject and template ' ...
        'fiducial triangles differ in shape. Check the selected points (S.scale = true allows ' ...
        'for a different head size).'], max(F.fid_error));
end

% Apply T to every template mesh
temp_brain = struct();
for k = 1:size(files, 1)
    m = meshes.(files{k, 1});
    m.vertices = apply_T(T, m.vertices);
    temp_brain.(files{k, 1}) = m;
end

if S.plot
    fid_reg = apply_T(T, fid_template);
    figure('Color', 'w'); hold on
    patch('Vertices', S.head.vertices, 'Faces', S.head.faces, ...
        'FaceColor', 'none', 'EdgeColor', 'k', 'EdgeAlpha', 0.1);
    patch('Vertices', temp_brain.brain.vertices, 'Faces', temp_brain.brain.faces, ...
        'FaceColor', 'none', 'EdgeColor', 'b', 'EdgeAlpha', 0.3);
    patch('Vertices', temp_brain.iskull.vertices, 'Faces', temp_brain.iskull.faces, ...
        'FaceColor', 'none', 'EdgeColor', 'r', 'EdgeAlpha', 0.3);
    patch('Vertices', temp_brain.scalp.vertices, 'Faces', temp_brain.scalp.faces, ...
        'FaceColor', 'none', 'EdgeColor', 'g', 'EdgeAlpha', 0.2);
    plot3(fid_subject(:, 1), fid_subject(:, 2), fid_subject(:, 3), 'o', ...
        'MarkerFaceColor', 'c', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
    plot3(fid_reg(:, 1), fid_reg(:, 2), fid_reg(:, 3), 's', ...
        'MarkerFaceColor', 'm', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
    text(fid_subject(:, 1), fid_subject(:, 2), fid_subject(:, 3), {'  NAS', '  LPA', '  RPA'});
    legend({'head', 'cortex', 'inner skull', 'scalp', 'subject fiducials', 'template fiducials'}, ...
        'Location', 'northeastoutside');
    title(sprintf('Brain registration: %.0f%% of cortex inside the head', 100 * F.frac_inside));
    axis equal; axis off; view(3);
    hold off
end
end


function F = fit_fiducials(fid_subject, fid_template, do_scale, brain_vertices, head)
% Units: a power of ten between subject and template (mm)
sf = 10 ^ round(log10(sqrt(tri_area(fid_subject) / tri_area(fid_template))));
s = sf;
if do_scale
    s = sqrt(tri_area(fid_subject) / tri_area(fid_template));
end
M0 = diag([s s s 1]);
fid_scaled = apply_T(M0, fid_template);

% Rigid fit: spm_eeg_inv_rigidreg takes points as COLUMNS (3 x n)
M1 = spm_eeg_inv_rigidreg(fid_subject', fid_scaled');
F.T = M1 * M0;
F.scale = s / sf;
F.fid_error = vecnorm(apply_T(F.T, fid_template) - fid_subject, 2, 2);

% Fraction of cortex inside the head (every 10th vertex)
pts = apply_T(F.T, brain_vertices(1:10:end, :));
F.frac_inside = mean(inside_mesh(pts, head.vertices, head.faces));
end


function a = tri_area(p)
a = norm(cross(p(1, :) - p(3, :), p(2, :) - p(3, :))) / 2;
end


function v = apply_T(T, v)
v = (T * [v, ones(size(v, 1), 1)]')';
v = v(:, 1:3);
end


function in = inside_mesh(pts, V, F)
% Winding number from the solid angle of every triangle (van Oosterom &
% Strackee 1983): about 4*pi inside a closed surface, about 0 outside.
% Orientation does not matter; a surface open far from the points (a torso
% cut at the waist) still gives the right answer near the head.
A = V(F(:, 1), :); B = V(F(:, 2), :); C = V(F(:, 3), :);
in = false(size(pts, 1), 1);
for i = 1:size(pts, 1)
    a = A - pts(i, :); b = B - pts(i, :); c = C - pts(i, :);
    la = vecnorm(a, 2, 2); lb = vecnorm(b, 2, 2); lc = vecnorm(c, 2, 2);
    num = dot(a, cross(b, c, 2), 2);
    den = la .* lb .* lc + dot(a, b, 2) .* lc + dot(a, c, 2) .* lb + dot(b, c, 2) .* la;
    in(i) = abs(sum(2 * atan2(num, den))) > 2 * pi;
end
end
