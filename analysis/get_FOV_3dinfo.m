function [output] = get_FOV_3dinfo( infile )
% usage: [output] = get_FOV_3dinfo( infile )
%
% input: infile = string, path to data file
%
% output: output = structure with fields:
%   - FOV = matrix, field of view, in mm
%   - vox_size = vector, voxel size, in mm
%   - n_vox = vector, number of voxels
%   - dim_labels = cell array of strings, dimension labels
%
% mps 20201215

%% load 3dinfo

[~, result] = system(['3dinfo ' infile]);

dim_labels = {'RL';'AP';'IS'};

key_pat_dim = {'R-to-L extent:', 'A-to-P extent:', 'I-to-S extent:'};
key_pat_FOV = {'] -to-' , '] -step-' };
key_pat_vox = {'mm [', 'voxels]'};

output.FOV = [];
output.vox_size = [];
output.n_vox = [];

for iDim = 1:numel(dim_labels)
    d_idx = strfind(result, key_pat_dim{iDim});
    
    get_line = result(d_idx+length(key_pat_dim{iDim}):d_idx+length(key_pat_dim{iDim})+100);
        
    %% first get FOV
    f1_idx = strfind(get_line, key_pat_FOV{1});
    
    if ~isempty(f1_idx)
        output.FOV(iDim,1) = str2num(...
            get_line( 1:(f1_idx(1)-3) ) ); % FOV 1 is before this str
    else
        error(['can''t find pattern' key_pat_FOV{1}]);
    end
    
    f2_idx = strfind(get_line, key_pat_FOV{2});
    
    if ~isempty(f2_idx)
        output.FOV(iDim,2) = str2num(...
            get_line( (f1_idx(1)+length(key_pat_FOV{1}) ):...
            (f2_idx(1)-3) ) ); % FOV 2 is between FOV 1 and here
    else
        error(['can''t find pattern' key_pat_FOV{2}]);
    end
    
    %% next get vox size and num
    v1_idx = strfind(get_line, key_pat_vox{1});
    
    if ~isempty(v1_idx)
        output.vox_size(iDim,1) = str2num(...
            get_line( (f2_idx(1)+length(key_pat_FOV{2}) ):...
            (v1_idx(1)-1) ) ); % voxel size is between FOV 2 and here
    else
        error(['can''t find pattern' key_pat_FOV{1}]);
    end
    
    v2_idx = strfind(get_line, key_pat_vox{2});
    
    if ~isempty(v2_idx)
        output.n_vox(iDim,1) = str2num(...
            get_line( (v1_idx(1)+length(key_pat_vox{1}) ):...
            (v2_idx(1)-1) ) ); % n voxels is between voxel size and here
    else
        error(['can''t find pattern' key_pat_FOV{1}]);
    end
    
end

%% out
output.dim_labels = dim_labels;
output.infile = infile;

end