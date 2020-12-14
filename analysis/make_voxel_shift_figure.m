function make_voxel_shift_figure( options )
% usage: make_voxel_shift_figure( options )
%
% input: options = structure with fields:
%   - ROI = string, 'vmPFC' or 'dmPFC'
%   - axial_slice_of_interest = integer
%   - img_scale = scalar, max and min for voxel shift maps
%   - topDir = string, path to top directory where analysis data live
%   - gitDir = string, path to git repo
%
% output: none (makes 4 figures)
%
% mps c. Sept. 2020

%% check matlab version
if datenum(version('-date')) < datenum('September 14, 2017')
    error(['This code relies on new functionality of niftiread.m, and will NOT function with versions '...
        'of Matlab older than 2017b'])
end

%% opts
if ~exist('options','var')
    options = [];
end
if ~isfield(options,'ROI')
    options.ROI = 'vmPFC';
end
if ~isfield(options,'axial_slice_of_interest')
    if strcmp(options.ROI,'vmPFC')
        options.axial_slice_of_interest = 22; % 22 for vmPFC, 46 for dmPFC
    elseif strcmp(options.ROI,'dmPFC')
        options.axial_slice_of_interest = 46; % 22 for vmPFC, 46 for dmPFC
    end
end
if ~isfield(options,'img_scale')
    options.img_scale = 10;
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

GE_dir = fullfile(subj_dir,'GE_qwarp');
SE_dir = fullfile(subj_dir,'SE_qwarp');
B0_dir = fullfile(subj_dir,'fugue');
raw_dir = fullfile(subj_dir,'raw');

EPI_file = fullfile(raw_dir,'PRF1_AP.nii.gz');
EPI_info = get_FOV_3dinfo(EPI_file);

ROI_file = fullfile(top_dir,[options.ROI '_ROI_mask.nii.gz']);
ROI_data = double(niftiread(ROI_file));

total_readout = 0.0416; % sec
% to convert ph_rad_per_sec to voxel shift -- multiply by total readout
% then divide by 2 pi

qwarp_file = 'blip_warp_For_WARP.nii.gz';
B0_file = 'ph_rad_per_sec_medfilt.nii.gz';

all_files.GE = fullfile(GE_dir, qwarp_file);
all_files.SE = fullfile(SE_dir, qwarp_file);
all_files.B0 = fullfile(B0_dir, B0_file);

PE_dimension = 2;

axial_slice_of_interest = options.axial_slice_of_interest;

sagittal_slice_of_interest = 72;

% dimensions are different for 3dQwarp and fugue voxel shift files,
% use the size of the EPI matrix to decide how big the images should be
% and center everything based on the position of the FOV (below)
max_dim = EPI_info.n_vox';

%% add path
addpath(options.gitDir);

%% first make sagittal GE EPI image -- do this in AFNI & save manually

%% get data & make voxel shift image
data_types = {'GE','B0','SE'};
add_title = {' oppPE',' FM',' oppPE'};
colors = {'r','g','b'};

max_line = 0;
min_line = 0;

for iImg = 1:numel(data_types)
%% load data
    clear load_data
    if ~exist(all_files.(data_types{iImg}),'file')
            brik_file = [all_files.(data_types{iImg})(1:end-7) '+orig.BRIK'];
            if ~exist(brik_file,'file')
                error(['Can''t find ' brik_file ' - have you run do_all_scans.sh?']);
            else
                [~, result] = system(['3dcopy ' brik_file ' ' ...
                    all_files.(data_types{iImg})]);
                load_data = double(niftiread(all_files.(data_types{iImg})));
            end
    else
        load_data = double(niftiread(all_files.(data_types{iImg})));
    end
    
    map_info = get_FOV_3dinfo(all_files.(data_types{iImg}));
    
    if ~strcmp(data_types{iImg},'B0')
        load_data = squeeze(load_data(:,:,:,1,PE_dimension));
    else
        load_data = load_data .* total_readout / 2 / pi;
    end
        
    FOV_diff = round( ( abs(map_info.FOV - EPI_info.FOV) )./repmat(...
        map_info.vox_size,[1 2]) );
    
    img_data = load_data( FOV_diff(1,1)+1 : end-FOV_diff(1,2) , ...
        FOV_diff(2,1)+1 : end-FOV_diff(2,2) , ...
        FOV_diff(3,1)+1 : end-FOV_diff(3,2) );
    
    if sum( size(img_data) == max_dim) ~= numel(max_dim)
        error(['cropped image size is not ' num2str(max_dim)]);
    end
    
    line_data.(data_types{iImg}) = img_data(end-sagittal_slice_of_interest,:,...
        axial_slice_of_interest);
    max_line = max([max_line max(line_data.(data_types{iImg})) ]);
    min_line = min([min_line min(line_data.(data_types{iImg})) ]);
    
%% plot voxel shift
    figure
    hold on
    
    plot_me = img_data(end:-1:1,end:-1:1,axial_slice_of_interest)';
    % end:-1:1 to flip so left = left, transpose so anterior = up
    
    plot_me = (plot_me + options.img_scale); % make so values < -options.img_scale are negative (min)
    plot_me = plot_me ./ (2*options.img_scale); % make so options.img_scale is max
    
    plot_me = repmat(plot_me, [1 1 3]);
    
    ROI_slice = ROI_data(end:-1:1,end:-1:1,axial_slice_of_interest)';
    ROI_idx = find(ROI_slice(:) == 1);
    [ROI_x, ROI_y] = ind2sub(size(ROI_slice),ROI_idx);
    
    for iIdx = 1:numel(ROI_x)
        plot_me(ROI_x(iIdx), ROI_y(iIdx), 2) = 0; % make ROI magenta
    end
    
    imagesc(plot_me,[0 1]);
    cmap = repmat([0:1/255:1]',[1 3]);
    plot([sagittal_slice_of_interest sagittal_slice_of_interest], ...
        [1 size(img_data,2)],'linestyle','-','color',colors{iImg},...
        'linewidth',2)
    colormap(cmap);
    box off
    cb = colorbar;
    set(cb,'color','k')
    cb.Label.String = 'Voxel shift';
    title([data_types{iImg} add_title{iImg}],'color','k');
    set(gca,'XTick',[],'YTick',[],'fontsize',18);
    set(gcf,'color','w')
    axis image
    axis off
    cb.Ticks = 0:0.25:1;
    cb.TickLabels = -options.img_scale:(options.img_scale*2/4):...
        options.img_scale;

%% plot line data
    if iImg == 1
        h = figure;
        hold on
        
        % plot a patch for where the ROI is
        ROI_idx = find(ROI_y == sagittal_slice_of_interest);
        ROI_min = min(ROI_x(ROI_idx));
        ROI_max = max(ROI_x(ROI_idx));
        patch([min_line-4 max_line+4 max_line+4 min_line-4], ...
            [ROI_max ROI_max ROI_min ROI_min], [1 0.9 1],...
            'LineStyle','none');
    else
        figure(h)
    end

    plot(line_data.(data_types{iImg})(end:-1:1),...
        1:numel(line_data.(data_types{iImg})),'linestyle','-',...
        'color',colors{iImg},'linewidth',2)
    if iImg == numel(data_types)
        set(gcf,'color','w')
        set(gca,'YTick',[1 numel(line_data.(data_types{iImg}))],...
            'YTickLabel',{'Pos.','Ant.'},...
            'XTick',-options.img_scale*2:2:options.img_scale*2,...
            'fontsize',18,'XColor','k','YColor','k')
        xlabel('Voxel shift (A-P)','color','k')
        axis([min_line-2 max_line+2 1 max_dim(2)]);
    end

end

end