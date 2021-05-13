function make_uncorr_sub_corr_figure( options )
% usage: make_uncorr_sub_corr_figure( options )
%
% input: options = structure with fields:
%   - axial_slice_of_interest = integer
%   - img_scale = scalar, max and min for image
%   - topDir = string, path to top directory where analysis data live
%
% output: none (makes 5 figures)
%
% mps 2021.02.02

%% check matlab version
if datenum(version('-date')) < datenum('September 14, 2017')
    error(['This code relies on new functionality of niftiread.m, and will NOT function with versions '...
        'of Matlab older than 2017b'])
end

%% opts
if ~exist('options','var')
    options = [];
end
if ~isfield(options,'axial_slice_of_interest')
    options.axial_slice_of_interest = 20;
end
if ~isfield(options,'img_scale')
    options.img_scale = 500;
end
if ~isfield(options,'topDir')
    options.topDir = input('Path to data directory: ','s');
end

%% settings
top_dir = options.topDir;
subj_dir = fullfile(top_dir, 'P1010228','B');
analyses = {'GE_qwarp', 'GE_topup', 'fugue', 'SE_qwarp', 'SE_topup'};
analyses_names = {'GE qwarp', 'GE topup', 'fugue', 'SE qwarp', 'SE topup'};

axial_slice = options.axial_slice_of_interest;

figure

for iA = 1:numel(analyses)
    %% load data
    sub_file = fullfile(subj_dir,analyses{iA},'uncorr_min_corr.nii.gz');
    if ~exist(sub_file,'file')
        error(['Can''t find ' sub_file ' - you probably need to create it manually...']);
        % e.g., 3dcalc -prefix filename.nii.gz -a 'file1.nii.gz' -b
        % 'file2.nii.gz' -expr 'a - b'
    else
        load_data = double(niftiread(sub_file));
    end
    
    %% plot
    subplot(2,3,iA+1); hold on
    
    
    plot_me = squeeze(load_data(end:-1:1, end:-1:1, axial_slice))';
    % end:-1:1 to flip so left = left, transpose so anterior = up
    
    plot_me = (plot_me + options.img_scale); % make so values < -options.img_scale are negative (min)
    plot_me = plot_me ./ (2*options.img_scale); % make so options.img_scale is max
    
    plot_me = repmat(plot_me, [1 1 3]);
    
    imagesc(plot_me,[0 1]);
    cmap = repmat([0:1/255:1]',[1 3]);
    
    colormap(cmap);
    box off
    set(gca,'XTick',[],'YTick',[],'fontsize',18);
    set(gcf,'color','w')
    axis image
    axis off
    title(analyses_names{iA})
end

end