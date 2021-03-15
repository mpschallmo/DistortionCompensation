function output = get_voxel_shift_data( options )
% usage: output = get_voxel_shift_data( options )
%
% mps 20210127
%% opt
if ~exist('options','var')
    options = [];
end
if ~isfield(options,'region')
    options.region = 'whole_brain';
    % valid options are: 'whole_brain', 'vmPFC', 'dmPFC', 'posterior'
end
if ~isfield(options,'subjDirs')
    options.subjDirs = {'P6003691'
        'P1010228'
        'P6001501'
        'P5104604'
        'P6004604'
        'P6004202'
        'P1010299'
        'P6004687'
        'P1007451'
        'P2104777'
        'P6010671'
        'P4100631'
        'P2110465'
        'P1010422'
        'P6010465'
        'P6004777'
        'P1010407'
        'P6010731'
        'P6010363'
        'P1006397'
        'P4110363'
        'P4104604'
        'P3110692'
        'P6010932'
        'P3102476'
        'P3111176'
        'P1011139'
        'P6004002'
        'P1011033'
        'P1010859'
        'P1011399'};
end
if ~isfield(options,'scanSubDirs')
    options.scanSubDirs = {'Z'
        'B'
        'B'
        'B'
        'Z'
        'Z'
        'B'
        'B'
        'B'
        'B'
        'B'
        'B'
        'B'
        'B'
        'Z'
        'Z'
        'B'
        'Z'
        'Z'
        'B'
        'B'
        'B'
        'B'
        'B'
        'B'
        'B'
        'B'
        'Z'
        'B'
        'B'
        'B'};
end

if ~isfield(options,'which_analysis')
    options.which_analysis = 'separate_GE'; % separate_GE = main analysis,
    % single_GE = analysis #2, SBRef = analysis #3
end
if ~isfield(options,'topDir')
    options.topDir = input('Path to data directory: ','s');
end
if ~contains(options.topDir, options.which_analysis)
    options.topDir = fullfile(options.topDir, options.which_analysis);
    warning(['Assuming you want to include ' options.which_analysis ...
        ' at the end of the data directory path: ' options.topDir]);
end
if ~isfield(options,'gitDir')
    options.gitDir = input('Path to git directory: ','s');
end
if ~contains(options.gitDir(end-8:end), 'analysis')
    options.gitDir = fullfile(options.gitDir, 'analysis');
    
    if ~exist(options.gitDir, 'dir')
        error(['Can''t find analysis directory in the specified location: '...
            options.gitDir]);
    end
end
addpath(genpath(options.gitDir));

if ~isfield(options,'overwrite_saved')
    options.overwrite_saved = 0; % 0 = no, 1 = yes
end
if ~isfield(options,'displayFigs')
    options.displayFigs = 1; % 0 = no, 1 = yes
end

output = [];

%% analysis
analyses = {'GE_qwarp', 'fugue', 'SE_qwarp'}; % only do qwarp, not topup
mat_file = fullfile(options.topDir, 'voxel_shift_data.mat');

qwarp_file = 'blip_warp_For_WARP.nii.gz';
B0_file = 'ph_rad_per_sec_medfilt.nii.gz';
T1_mask = '3T_anat_uni_al_EPI_mask.nii.gz';
EPI_files = {'PRF1_AP.nii.gz','CSS_task3_AP.nii.gz','COP_task1_AP.nii.gz'};
GU_file = 'fullWarp_rel.nii.gz';

PE_dimension = 2;
total_readout = 0.0416; % sec
% to convert ph_rad_per_sec to voxel shift -- multiply by total readout
% then divide by 2 pi

if ~strcmp(options.region, 'whole_brain') % if using ROI mask
    ROI_file = fullfile(options.topDir,[options.region '_ROI_mask.nii.gz']);
    ROI_data = double(niftiread(ROI_file));
    add_GU = 0;
else
    add_GU = 1;
end

if exist(mat_file,'file') % file exists
    load_data = load(mat_file);
    output = load_data.output;
    if ~options.overwrite_saved && isfield(output, options.region)
        overwrite = 0; % there is a structure for this region in the file, and not trying to overwrite, so use it
        
        % check if # of subjects in saved data file matches options
        if numel(options.subjDirs) ~= numel(output.(options.region...
                ).(analyses{1}).subj_data.mean)
            error(['Number of subjects in options.subjDirs = ' num2str(numel(options.subjDirs))...
                ' but number in saved data file = ' num2str(numel(output.(options.region...
                ).(analyses{1}).subj_data.mean)) ' ... perhaps you want to overwrite?']);
        end
        
    else
        overwrite = 1;
    end
else
    overwrite = 1; % no file saved, or we said we want to overwrite
end

if overwrite
    h_wait = waitbar(0, 'loading data, please wait...');
    for iAnalysis = 1:numel(analyses)
        for iSubj = 1:numel(options.subjDirs)
            waitbar((iSubj + (iAnalysis-1)*numel(options.subjDirs) )./...
                ((numel(analyses)+add_GU)*numel(options.subjDirs)), h_wait);
            if strcmp(analyses{iAnalysis},'fugue')
                use_map = B0_file;
            else
                use_map = qwarp_file;
            end
            subj_vox_map_file = fullfile(options.topDir, options.subjDirs{iSubj},...
                options.scanSubDirs{iSubj}, analyses{iAnalysis}, use_map);
            if ~exist(subj_vox_map_file,'file')
                convert_BRIK(subj_vox_map_file) % function @ bottom
            end
            
            epi_file = fullfile(options.topDir, options.subjDirs{iSubj},...
                options.scanSubDirs{iSubj}, 'raw', EPI_files{iAnalysis});
            T1_mask_file = fullfile(options.topDir, options.subjDirs{iSubj},...
                options.scanSubDirs{iSubj}, analyses{iAnalysis}, T1_mask);
            
            clear load_data
            load_data = double(niftiread(subj_vox_map_file));
            map_info = get_FOV_3dinfo(subj_vox_map_file);
            
            EPI_info = get_FOV_3dinfo(epi_file);
            
            if ~strcmp(analyses{iAnalysis},'fugue')
                load_data = squeeze(load_data(:,:,:,1,PE_dimension));
            else
                load_data = load_data .* total_readout / 2 / pi;
            end
            
            FOV_diff = (map_info.FOV - EPI_info.FOV)./repmat(...
                map_info.vox_size,[1 2]);
            
            if sum(round(FOV_diff(:,1) - FOV_diff(:,2))) == 0 % this is a shift, not a size difference
                if sum(round(FOV_diff) > 1)
                    error('more than 1 voxel shift between EPI and FM...')
                else
                    img_data = load_data;
                end
            else
                FOV_diff = round( ( abs(map_info.FOV - EPI_info.FOV) )./repmat(...
                    map_info.vox_size,[1 2]) );
                img_data = load_data( FOV_diff(1,1)+1 : end-FOV_diff(1,2) , ...
                    FOV_diff(2,1)+1 : end-FOV_diff(2,2) , ...
                    FOV_diff(3,1)+1 : end-FOV_diff(3,2) ); % crop to original size, 130 x 130 x 85
            end
            
            hist_data = double(niftiread(T1_mask_file)); % start with T1 mask
            
            if ~strcmp(options.region, 'whole_brain') % if using ROI mask
                hist_data = hist_data .* ROI_data; % mask with ROI
            end
            
            hist_data( hist_data == 0 ) = NaN; % remove data outside mask
            
            hist_data = abs( img_data .* hist_data ); % mask the voxel map data, take absolute value
            
            output.(options.region).(analyses{iAnalysis}).subj_data.mean(iSubj) = ...
                nanmean(hist_data(:)); % mean across voxels, ignore NaNs
            output.(options.region).(analyses{iAnalysis}).subj_data.median(iSubj) = ...
                nanmedian(hist_data(:)); % median across voxels, ignore NaNs
            output.(options.region).(analyses{iAnalysis}).subj_data.sd(iSubj) = ...
                nanstd(hist_data(:)); % sd across voxels, ignore NaNs
            output.(options.region).(analyses{iAnalysis}).subj_data.max(iSubj) = ...
                nanmax(hist_data(:)); % max across voxels, ignore NaNs
        end
        output.(options.region).(analyses{iAnalysis}).group_mean_sd.mean = ...
            [mean(output.(options.region).(analyses{iAnalysis}).subj_data.mean) ...
            std(output.(options.region).(analyses{iAnalysis}).subj_data.mean)];
        
        output.(options.region).(analyses{iAnalysis}).group_mean_sd.median = ...
            [mean(output.(options.region).(analyses{iAnalysis}).subj_data.median) ...
            std(output.(options.region).(analyses{iAnalysis}).subj_data.median)];
        
        output.(options.region).(analyses{iAnalysis}).group_mean_sd.sd = ...
            [mean(output.(options.region).(analyses{iAnalysis}).subj_data.sd) ...
            std(output.(options.region).(analyses{iAnalysis}).subj_data.sd)];
        
        output.(options.region).(analyses{iAnalysis}).group_mean_sd.max = ...
            [mean(output.(options.region).(analyses{iAnalysis}).subj_data.max) ...
            std(output.(options.region).(analyses{iAnalysis}).subj_data.max)];
        
    end
    
    %% include gradunwarp
    if strcmp(options.region , 'whole_brain')
        for iSubj = 1:numel(options.subjDirs)
            waitbar((iSubj + numel(analyses)*numel(options.subjDirs) )./...
                ((numel(analyses)+add_GU)*numel(options.subjDirs)), h_wait);
            GU_data_file = fullfile(options.topDir, options.subjDirs{iSubj},...
                options.scanSubDirs{iSubj}, 'raw', GU_file);
            
            clear img_data
            img_data = double(niftiread(GU_data_file));
            
            hist_data = double(niftiread(T1_mask_file)); % start with T1 mask
            
            hist_data( hist_data == 0 ) = NaN; % remove data outside mask
            
            for iDim = 1:size(img_data,4)
                hist_data = abs( squeeze(img_data(:,:,:,iDim)) .* hist_data ); % mask the voxel map data, take absolute value
                
                output.(options.region).gradunwarp(iDim).subj_data.mean(iSubj) = ...
                    nanmean(hist_data(:)); % mean across voxels, ignore NaNs
                output.(options.region).gradunwarp(iDim).subj_data.median(iSubj) = ...
                    nanmedian(hist_data(:)); % median across voxels, ignore NaNs
                output.(options.region).gradunwarp(iDim).subj_data.sd(iSubj) = ...
                    nanstd(hist_data(:)); % sd across voxels, ignore NaNs
                output.(options.region).gradunwarp(iDim).subj_data.max(iSubj) = ...
                    nanmax(hist_data(:)); % max across voxels, ignore NaNs
            end
        end
        for iDim = 1:size(img_data,4)
            output.(options.region).gradunwarp(iDim).group_mean_sd.mean = ...
                [mean(output.(options.region).gradunwarp(iDim).subj_data.mean) ...
                std(output.(options.region).gradunwarp(iDim).subj_data.mean)];
            
            output.(options.region).gradunwarp(iDim).group_mean_sd.median = ...
                [mean(output.(options.region).gradunwarp(iDim).subj_data.median) ...
                std(output.(options.region).gradunwarp(iDim).subj_data.median)];
            
            output.(options.region).gradunwarp(iDim).group_mean_sd.sd = ...
                [mean(output.(options.region).gradunwarp(iDim).subj_data.sd) ...
                std(output.(options.region).gradunwarp(iDim).subj_data.sd)];
            
            output.(options.region).gradunwarp(iDim).group_mean_sd.max = ...
                [mean(output.(options.region).gradunwarp(iDim).subj_data.max) ...
                std(output.(options.region).gradunwarp(iDim).subj_data.max)];
        end
    end
    %% save data
    save(mat_file,'output');
    close(h_wait);
end

%% plot
if options.displayFigs
    
    h = figure;
    use_symb = {'s','^','s'};
    use_color = {'r',[0 0.75 0],'b'};
    use_fill = {'r',[0 0.75 0],'b'};
    use_mk_size = 10;
    use_lw = 1;
    use_x_labels = {'GE','B0','SE'};
    use_y_labels = {'mean','median','S.D.',...
        'max.'};
    use_title = [options.region ' voxel shift'];
    
    for iAnalysis = 1:numel(analyses)
        
        plot_me  = [output.(options.region).(analyses{iAnalysis}).group_mean_sd.mean ;
            output.(options.region).(analyses{iAnalysis}).group_mean_sd.median ;
            output.(options.region).(analyses{iAnalysis}).group_mean_sd.sd ;
            output.(options.region).(analyses{iAnalysis}).group_mean_sd.max];
        figure(h)
        for iP = 1:size(plot_me,1)
            subplot(4,1,iP)
            hold on
            errorbar(iAnalysis, plot_me(iP,1), plot_me(iP,2),...
                'Marker', use_symb{iAnalysis}, 'color', use_color{iAnalysis},...
                'MarkerFaceColor', use_fill{iAnalysis}, ...
                'MarkerSize', use_mk_size,'linewidth',use_lw)
            box off
            
            if iAnalysis == 3
                set(gca,'XTick',1:3,'XTickLabel',use_x_labels,...
                    'XColor','k','YColor','k','fontsize',18)
                ax = axis;
                if ax(3) < 0
                    ax(3) = -1;
                end
                axis([0.5 3.5 ax(3) ax(4)]);
                set(gcf,'color','w','POS',[5 300 330 665],'renderer','painters')
                ylabel(use_y_labels{iP},'color','k');
                if iP == 1
                   title(use_title,'color','k') 
                end
            end
        end
    end
end

%% out
output.options = options;

%% convert to .nii.gz function
    function result = convert_BRIK( filename )
        brik_file = [filename(1:end-7) '+orig.BRIK'];
        if ~exist(brik_file,'file')
            error(['Can''t find ' brik_file ' - have you run do_all_scans.sh?']);
        else
            [~, result] = system(['3dcopy ' brik_file ' ' filename]);
        end
    end
%%
end