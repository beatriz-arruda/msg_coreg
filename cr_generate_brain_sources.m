% cr_generate_brain_sources - Brain source model from the registered head meshes
%
% Builds the brain source model used for forward modelling of the cortical
% response, from the SPM template meshes that cr_register_brain has placed
% in sensor space. Two layouts are available:
%
%   'grid'    a regular 3D grid inside the inner skull (default). The inner
%             skull mesh is closed, so every point can be tested for being
%             inside it. This is the usual source space for MEG
%             beamforming.
%   'cortex'  the vertices of the cortical surface mesh, optionally
%             reduced to a target number, with their surface normals. The
%             template cortex is not a closed surface, so it is only used
%             for its vertices, never as a volume.
%
% USAGE:
%   src = cr_generate_brain_sources(S)
%
% INPUT:
%   S               - Structure with the following fields:
%
%   Required:
%     S.iskull        - inner skull mesh (.vertices, .faces), for 'grid'
%     S.brain         - cortical mesh (.vertices, .faces), for 'cortex'
%
%   Optional:
%     S.mode          - 'grid' (default) or 'cortex'
%     S.resolution    - grid spacing in mm (default: 5)
%     S.inwardshift   - distance kept from the inner skull in mm
%                       (default: 5); keeps sources away from the boundary,
%                       where BEM accuracy is lowest
%     S.n_vertices    - target number of cortical vertices for 'cortex'
%                       (default: all)
%     S.unit          - units of the meshes (default: 'mm')
%
% OUTPUT:
%   src             - Source model struct:
%                       .pos     [N x 3] source positions
%                       .inside  logical [N x 1] (all true)
%                       .unit    units (as S.unit)
%                       .mode    'grid' or 'cortex'
%                       .nrm     [N x 3] outward normals ('cortex' only)
%                       .resolution  grid spacing ('grid' only)
%
% SAVE IT IN THE GEOMETRY FILE AS sources_brain, next to mesh_brain,
% mesh_iskull, mesh_oskull and mesh_scalp. msg_fwd then detects the brain
% and computes brain lead fields (three-shell BEM and Biot-Savart) along
% with the cord.
%
% NOTES:
%   - A 5 mm grid in an adult inner skull gives roughly 8,000-10,000
%     sources; source reconstruction time grows with this number.
%   - The grid uses ft_prepare_sourcemodel with the inner skull as a
%     single-shell boundary.
%
% EXAMPLE:
%   all_meshes = cr_check_registration(S);           % with S.brain = true
%   B = [];
%   B.iskull     = all_meshes.iskull;
%   B.resolution = 5;
%   sources_brain = cr_generate_brain_sources(B);
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

function src = cr_generate_brain_sources(S)

if ~isfield(S, 'mode'),        S.mode = 'grid';   end
if ~isfield(S, 'resolution'),  S.resolution = 5;  end
if ~isfield(S, 'inwardshift'), S.inwardshift = 5; end
if ~isfield(S, 'n_vertices'),  S.n_vertices = []; end
if ~isfield(S, 'unit'),        S.unit = 'mm';     end

switch S.mode
    case 'grid'
        if ~isfield(S, 'iskull')
            error('cr_generate_brain_sources: S.iskull is required for the grid layout.');
        end
        bnd = struct('pos', double(S.iskull.vertices), 'tri', double(S.iskull.faces), ...
            'unit', S.unit);
        vol = struct('bnd', bnd, 'type', 'singleshell', 'unit', S.unit);
        cfg = [];
        cfg.headmodel   = vol;
        cfg.resolution  = S.resolution;
        cfg.unit        = S.unit;
        cfg.inwardshift = S.inwardshift;
        sm = ft_prepare_sourcemodel(cfg);
        src.pos = sm.pos(sm.inside, :);
        src.resolution = S.resolution;

    case 'cortex'
        if ~isfield(S, 'brain')
            error('cr_generate_brain_sources: S.brain is required for the cortex layout.');
        end
        mesh = struct('vertices', double(S.brain.vertices), 'faces', double(S.brain.faces));
        if ~isempty(S.n_vertices) && S.n_vertices < size(mesh.vertices, 1)
            mesh = reducepatch(mesh, S.n_vertices / size(mesh.vertices, 1));
        end
        g = gifti(struct('vertices', mesh.vertices, 'faces', mesh.faces));
        nrm = spm_mesh_normals(g, true);   % vertex normals (first output)
        src.pos = mesh.vertices;
        src.nrm = nrm;

    otherwise
        error('cr_generate_brain_sources: unknown mode ''%s''.', S.mode);
end

src.inside = true(size(src.pos, 1), 1);
src.unit   = S.unit;
src.mode   = S.mode;
fprintf('cr_generate_brain_sources: %d sources (%s)\n', size(src.pos, 1), S.mode);
end
