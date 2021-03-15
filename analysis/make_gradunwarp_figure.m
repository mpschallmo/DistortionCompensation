function make_gradunwarp_figure( options )
% usage: make_gradunwarp_figure( options )
%
% input: options = structure with fields:
%   - axial_slice_of_interest = integer
%   - coronal_slice_of_interest = integer
%   - sagittal_slice_of_interest = integer
%   - img_scale = scalar, max and min for voxel shift maps
%   - topDir = string, path to top directory where analysis data live
%   - gitDir = string, path to git repo
%
% output: none (makes 3x3 figures)
%
% mps c. Feb 2021

%% check matlab version
if datenum(version('-date')) < datenum('September 14, 2017')
    error(['This code relies on new functionality of niftiread.m, and will NOT function with versions '...
        'of Matlab older than 2017b'])
end

%% opts
if ~exist('options','var')
    options = [];
end
if ~isfield(options,'sagittal_slice_of_interest')
    options.sagittal_slice_of_interest = 58;
end
if ~isfield(options,'coronal_slice_of_interest')
    options.coronal_slice_of_interest = 60;
end
if ~isfield(options,'axial_slice_of_interest')
    options.axial_slice_of_interest = 43;
end
if ~isfield(options,'img_scale')
    options.img_scale = 4;
end
if ~isfield(options,'show_slices')
    options.show_slices = 0;
end
if ~isfield(options,'show_brain_mask')
    options.show_brain_mask = 1;
end
if ~isfield(options,'topDir')
    options.topDir = input('Path to data directory: ','s');
end
if ~isfield(options,'gitDir')
    options.gitDir = input('Path to git directory: ','s');
end

%% settings
top_dir = options.topDir;
subj_dir = fullfile(top_dir, 'P1006397','B');

raw_dir = fullfile(subj_dir,'raw');
uncorr_dir = fullfile(subj_dir,'uncorr');

GU_file = fullfile(raw_dir,'fullWarp_rel.nii.gz');
T1_mask_file = fullfile(uncorr_dir,'3T_anat_uni_al_EPI_mask.nii.gz');

sagittal_slice = options.sagittal_slice_of_interest;
coronal_slice = options.coronal_slice_of_interest;
axial_slice = options.axial_slice_of_interest;

%% add path
addpath(options.gitDir);

%% load data
if ~exist(GU_file,'file')
    error(['Can''t find ' GU_file ' - have you run do_all_scans.sh?']);
else
    load_data = double(niftiread(GU_file));
end

if ~exist(T1_mask_file,'file')
    error(['Can''t find ' T1_mask_file ' - have you run do_all_scans.sh?']);
else
    T1_data = double(niftiread(T1_mask_file));
end

%% plot voxel shift
colors = {'r','g','b'}; % sagittal, coronal, axial
titles = {'x','y','z';'sagittal','coronal','axial'};

figure

for iShiftDim = 1:size(load_data,4) % shifts in x, y, z stored in 4th dim
    for iView = 1:3 % sagittal, coronal, axial views
        subplot(3, 3, iView + 3*(iShiftDim-1))
        hold on
        
        if iView == 1
            plot_me = squeeze(load_data(sagittal_slice, 1:end, 1:end, iShiftDim))';
        elseif iView == 2
            plot_me = squeeze(load_data(end:-1:1, coronal_slice, 1:end, iShiftDim))';
        elseif iView == 3
            plot_me = squeeze(load_data(end:-1:1, end:-1:1, axial_slice, iShiftDim))';
            % end:-1:1 to flip so left = left, transpose so anterior = up
        end
        
        plot_me = (plot_me + options.img_scale); % make so values < -options.img_scale are negative (min)
        plot_me = plot_me ./ (2*options.img_scale); % make so options.img_scale is max
        
        plot_me = repmat(plot_me, [1 1 3]);
        
        if options.show_brain_mask
            if iView == 1
                T1_slice = squeeze(T1_data(sagittal_slice, 1:end, 1:end))';
            elseif iView == 2
                T1_slice = squeeze(T1_data(end:-1:1, coronal_slice, 1:end))';
            elseif iView == 3
                T1_slice = squeeze(T1_data(end:-1:1,end:-1:1,axial_slice))';
                % end:-1:1 to flip so left = left, transpose so anterior = up
            end
            T1_idx = find(T1_slice(:) == 1);
            [T1_x, T1_y] = ind2sub(size(T1_slice),T1_idx);

            for iIdx = 1:numel(T1_x)
                plot_me(T1_x(iIdx), T1_y(iIdx), 2:3) = 0; % make ROI red
            end
        end
        
        imagesc(plot_me,[0 1]);
        cmap = repmat([0:1/255:1]',[1 3]);
        
        if options.show_slices
            if iView == 1
                plot([coronal_slice coronal_slice], [1 size(load_data,3)], ...
                    'linestyle','--','color',colors{2},...
                    'linewidth',2)
                plot([1 size(load_data,2)], [axial_slice axial_slice], ...
                    'linestyle','--','color',colors{3},...
                    'linewidth',2)
            elseif iView == 2
                plot([sagittal_slice sagittal_slice], [1 size(load_data,3)], ...
                    'linestyle','--','color',colors{1},...
                    'linewidth',2)
                plot([1 size(load_data,1)], [axial_slice axial_slice], ...
                    'linestyle','--','color',colors{3},...
                    'linewidth',2)
            elseif iView == 3
                plot([sagittal_slice sagittal_slice], [1 size(load_data,2)], ...
                    'linestyle','--','color',colors{1},...
                    'linewidth',2)
                plot([1 size(load_data,1)], [coronal_slice coronal_slice], ...
                    'linestyle','--','color',colors{2},...
                    'linewidth',2)
            end
        end
        
        colormap(cmap);
        box off
        cb = colorbar;
        set(cb,'color','k')
        cb.Label.String = 'Voxel shift';
        title([titles{1,iShiftDim} ', ' titles{2,iView}],'color','k');
        set(gca,'XTick',[],'YTick',[],'fontsize',18);
        set(gcf,'color','w')
        axis image
        axis off
        cb.Ticks = 0:0.25:1;
        cb.TickLabels = -options.img_scale:(options.img_scale*2/4):...
            options.img_scale;
    end
end

end