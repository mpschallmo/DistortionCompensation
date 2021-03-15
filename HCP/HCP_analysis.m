function results = HCP_analysis( options )
% usage: results = HCP_analysis( options )
%
% input: options = structure with fields:
%   - subjDirs = cell array of strings, with subject IDs
%   - topDir = string, path to top directory where analysis data live
%   - displayFigs = binary, 1 = show figures, 0 = hide
%   - which_error_bars = string, options are Morey (within-subjects error
%       bars, as described by Morey (2008), sem (standard error of the mean)
%   - connect_not_sig = binary, 1 = show gray lines connecting
%       non-significant differences, 0 = no
%   - overwrite_data_file = binary, whether or not to overwrite mat file
%       in which analysis data are stored, 1 = yes, 0 = no
%   - show_anovas = binary, show figures for anovas, 1 = yes, 0 = no
%   - gitDir = string, path to git repo
%
% output: results = structure with fields:
%   - subj_number = numeric, subject ID #
%   - anova = structure with results from anovas comparing results between
%       conditions
%   - ttest = structure with results from paired t-tests between conditions
%   - save_FDR_p = cell array of matricies of FDR corrected p-values
%   - dice = structure with Dice coefficent data
%   - use_dice = matrix with Dice coefficient data
%   - use_minfo = matrix with mutual info data
%   - options = structure with input options
%   - date_run = string, time stamp for when the analysis was run
%
% mps c. July 2019

%% check matlab version
if datenum(version('-date')) < datenum('September 14, 2017')
    error(['This code relies on new functionality of niftiread.m, and will NOT function with versions '...
        'of Matlab older than 2017b'])
end

%% opts
if ~exist('options','var')
    options = [];
end
if ~isfield(options,'subjDirs')
    options.subjDirs = {'100610'
        '102311'
        '102816'
        '104416'
        '105923'
        '108323'
        '109123'
        '111312'
        '111514'
        '114823'
        '115017'
        '115825'
        '116726'
        '118225'
        '125525'
        '126426'
        '126931'
        '128935'
        '130114'
        '130518'};
end

if ~isfield(options,'topDir')
    options.topDir = input('Path to data directory: ','s');
end
if ~contains(options.topDir, 'HCP')
    options.topDir = fullfile(options.topDir, 'HCP');
    warning(['Assuming you want to include HCP' ...
        ' at the end of the data directory path: ' options.topDir]);
end
if ~isfield(options,'displayFigs')
    options.displayFigs = 1; % 1 = show, 0 = hide
end
if ~isfield(options,'which_error_bars')
    options.which_error_bars = 'Morey'; % Morey, sem
end
if ~isfield(options, 'connect_not_sig')
    options.connect_not_sig = 1; % gray lines connecting non-significant differences
end
if ~isfield(options, 'overwrite_data_file')% mps 20201102
    options.overwrite_data_file = 0;
end
if ~isfield(options, 'show_anovas')
    options.show_anovas = 0;
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

results = [];

%% get subject numbers
results.subj_number = [];
for iS = 1:numel(options.subjDirs)
    results.subj_number(iS) = str2num(options.subjDirs{iS});
end

%% load data
analyses = {'SE_qwarp','SE_topup','uncorr', '12param'};
analy_names = {'SE_qwarp', 'SE_topup', 'uncorr', 'param12'};
legend_names = {'SE 3dQwarp', 'SE topup', '6 param.', ...
    '12 param.'};

use_symb = {'s','^','o','o'};
use_color = {'b','b','k','k'};
use_fill = {'b','b','w',[0.5 0.5 0.5]};
quant_names = {'mask','mut_info'};
plot_q_names = {'Dice: whole-brain masks','Mutual info.: masked EPI & T1 data'};

nSubj = numel(options.subjDirs);
for iA = 1:numel(analyses)
    dice.(analy_names{iA}).mask = nan(nSubj,1);
end

h_wait = waitbar(0, 'loading FMSD data, please wait...');

data_file = fullfile(options.topDir, 'dist_comp_data.mat');

if exist(data_file,'file') % mps 20201102
    get_data = load(data_file); % this will preserve the other ROI files when overwriting
    data = get_data.data;
end

if options.overwrite_data_file || ~exist(data_file,'file')    
    % mps 20201102 make this a structure...
    data.dice = nan(nSubj,numel(analyses),numel(fieldnames(dice.(analy_names{1}))) );
    data.minfo = nan(nSubj,numel(analyses));    
    
    for iS = 1:nSubj
        subjDir = fullfile(options.topDir,options.subjDirs{iS},'proc');
        % n.b. assumes this directory structure, set by align_data.sh
        
        for iA = 1:numel(analyses)
            % get Dice data
            
            anat_mask = fullfile(subjDir,analyses{iA},['T1_' analyses{iA} ...
                '_mask.nii.gz']);
            if ~exist(anat_mask,'file')
                error(['Can''t find file: ' anat_mask]);
            end
            epi_mask = fullfile(subjDir,analyses{iA},['epi_' analyses{iA} ...
                '_mask.nii.gz']);
            if ~exist(epi_mask,'file')
                error(['Can''t find file: ' epi_mask]);
            end
            
            [~, result] = system(['3ddot -dodice ' anat_mask ' '...
                epi_mask]);
            
            num_idx = strfind(result,'0.');
            result = result(num_idx(end):end); % cut out warning message
            
            data.dice(iS,iA,1) = str2num(result);
            
            dice.(analy_names{iA}).mask(iS) = data.dice(iS,iA,1);
                        
            anat_3T_file = fullfile(subjDir,analyses{iA},['T1_' analyses{iA} ...
                '.nii.gz']);
            if ~exist(anat_3T_file,'file')
                error(['Can''t find file: ' anat_3T_file]);
            end
            epi_7T_file = fullfile(subjDir,analyses{iA},['epi_' analyses{iA} ...
                '.nii.gz']);
            if ~exist(epi_7T_file,'file')
                error(['Can''t find file: ' epi_7T_file]);
            end
            
            anat_3T_data = round(double(niftiread(anat_3T_file)));
            epi_7T_data = round(double(niftiread(epi_7T_file)));
            
            anat_3T_mask = round(double(niftiread(anat_mask)));
            epi_7T_mask = round(double(niftiread(epi_mask)));
            
            % mask mutual info, brain only
            anat_3T_data = anat_3T_data .* anat_3T_mask;
            epi_7T_data = epi_7T_data .* epi_7T_mask;
            
            data.minfo(iS,iA) = mutInfo(...
                anat_3T_data(:), epi_7T_data(:));
                
        end
        waitbar(iS/nSubj, h_wait);
    end
    
    save(data_file, 'data'); % mps 20201102
end
close(h_wait);

use_minfo = data.minfo;
use_dice = data.dice;

%% run some stats
if options.show_anovas
    showFig = 'on';
else
    showFig = 'off';
end

save_FDR_p = [];

cond_idx = [1 2]; % mask, mututal info

% first do an anova for each tissue type, all conditions
for use_idx = cond_idx
    if use_idx < 2
        all_data = use_dice(:,:,use_idx);
    elseif use_idx == 2
        all_data = use_minfo;
    else
        error(['unknown value of use_idx = ' num2str(use_idx) ])
    end
    
    all_subj = repmat([1:size(all_data,1)]',[1 size(all_data,2)]);
    all_cond = repmat([1:size(all_data,2)],[size(all_data,1) 1]);
    
    [panova, results.anova.(quant_names{use_idx}).all_con, stats] = anovan(all_data(:), ...
        {all_subj(:),  all_cond(:)}, 'random', 1, 'model', 'full', 'varnames', ...
        {'subj', 'condition'}, 'display', showFig);
    
    results.ttest.(quant_names{use_idx}) = cell(...
        (size(use_dice,2)-1), size(use_dice,2) );
    all_p = [];
    all_iC = [];
    save_FDR_p{use_idx} = nan((size(use_dice,2)-1), size(use_dice,2));
    % 1 tailed t-tests for all paired comparisons
    for iC1 = 1:(size(use_dice,2)-1)
        for iC2 = (iC1+1):size(use_dice,2)
            [h, p, ci, results.ttest.(quant_names{use_idx}){iC1,iC2}] = ...
                ttest(all_data(:, iC1), all_data(:, iC2));
            results.ttest.(quant_names{use_idx}){iC1,iC2}.p_uncorr = p;
            all_p = [all_p ; p];
            all_iC = [all_iC ; iC1 iC2];
            results.ttest.(quant_names{use_idx}){iC1,iC2}.C1_name = ...
                analy_names{iC1};
            results.ttest.(quant_names{use_idx}){iC1,iC2}.C2_name = ...
                analy_names{iC2};
        end
    end
    [sort_p, p_idx] = sort(all_p, 'ascend');
    FDR_p = sort_p.*[numel(sort_p):-1:1]';
    for iP = 1:numel(FDR_p)
        results.ttest.(quant_names{use_idx}){all_iC(p_idx(iP), 1), ...
            all_iC(p_idx(iP), 2)}.p_FDR = ...
            FDR_p(iP);
        save_FDR_p{use_idx}(all_iC(p_idx(iP), 1), ...
            all_iC(p_idx(iP), 2)) = FDR_p(iP);
    end
    
end

%% display figs
subplot_vals = {[0.10 0.71 0.88 0.25],...
    [0.10 0.38 0.88 0.25]}; % used to set size of plot regions

subplot_vals2 = {[0.09 0.57 0.39 0.38],...
    [0.59 0.57 0.39 0.38]};

use_x_labels = {'SE oppPE','Align only'};
use_x_tick = [1.5 3.5];
use_mk_size = 10;
use_lw = 1;
use_patch_pct_y = 0.05;
use_patch_pct_x = 0.01;

if options.displayFigs
    h1 = figure;
    h2 = figure;
    
    for iM = cond_idx
        % box plots
        figure(h1)
        subplot('position',subplot_vals{iM}); hold on
        plot([0 0],[0 0],'r-','linewidth',use_lw)
        plot([0 0],[0 0],'b-','linewidth',use_lw)
        plot([0 0],[0 0],'k--','linewidth',use_lw)
        
        if iM < 2
            plot_me = squeeze(use_dice(:,:,iM));
        else
            plot_me = use_minfo;
        end
        which_analysis = repmat([1:numel(analyses)], ...
            [size(use_dice,1) 1]);
        hb = boxplot(plot_me(:), which_analysis(:));
        hold on
        set(hb,'linewidth',use_lw)
        bee_bin_width = 0.075;
        bee_spread_width = 0.9;
        pause(0.2)
        
        for iA = 1:numel(analyses)
            figure(h1)
            if iM < 2
                mask_data = use_dice(:,iA,iM);
            else
                mask_data = use_minfo(:,iA);
            end
            hp = plotSpread({mask_data},'binWidth',...
                bee_bin_width*nanmean(mask_data),...
                'distributionColors',{[0.7 0.7 0.7]},...
                'xValues',iA,'spreadWidth',bee_spread_width);
            set(hp{1},'MarkerSize',14)
        end
        hold on
        
        box off
        title(plot_q_names{iM})
        set(gca,'XTick',use_x_tick,'XTickLabel',use_x_labels,...
            'fontsize',18,'XColor','k','YColor','k')
        set(gcf,'color','w')
        if iM < 2
            ylabel('Dice coeff.','color','k')
        else
            ylabel('Mutual info.','color','k')
        end
        ax = axis;
        axis([0.5 numel(analyses)+0.5 ax(3) ax(4)]);
        
        patch_x = use_patch_pct_x*(numel(analyses)+0.5 - 0.5);
        patch_y = use_patch_pct_y*(ax(4) - ax(3));
        patch([0.5-patch_x 0.5+patch_x 0.5+patch_x 0.5-patch_x],...
            [ax(3)-patch_y ax(3)-patch_y ax(3)+patch_y ax(3)+patch_y],...
            'w','EdgeColor','none');
        
        pause(0.2)
        
    end
    
    
    for iA = 1:size(use_dice,2)
        for iM = cond_idx
            % means + error bars
            
            if strcmp(lower(options.which_error_bars),'sem') % show SEM error bars
                if iM < 2
                    use_data = use_dice(:,:,iM);
                    plot_me = use_dice(:,iA,iM);
                else
                    use_data = use_minfo;
                    if iA == 1 % for plotting mut info 
                        plot_me = use_minfo(:,iA);
                    else
                        plot_me = use_minfo(:,iA);
                    end
                end
            else strcmp(lower(options.which_error_bars),'morey') % show morey within-subj error bars
                if iM < 2
                    use_data = use_dice(:,:,iM) - ...
                        repmat(nanmean(use_dice(:,:,iM),2),[1 size(use_dice,2)]) + ...
                        repmat(nanmean(nanmean(use_dice(:,:,iM),2),1), ...
                        [size(use_dice,1) size(use_dice,2)]); % subtract individual subject mean, add group mean
                    plot_me = use_data(:,iA);
                else
                    use_data = use_minfo - repmat(nanmean(use_minfo,2),...
                        [1 size(use_minfo,2)]) + repmat(nanmean(nanmean(...
                        use_minfo,2),1), size(use_minfo));
                    if iA == 1 % for plotting mut info
                        plot_me = use_data(:,iA);
                    else
                        plot_me = use_data(:,iA);
                    end
                end
            end
            
            figure(h2)
            subplot('position',subplot_vals2{iM}); hold on
            if iM == 2 && iA == 1
                for iS = 1:numel(use_symb)
                    plot(-1,-1,'Marker',use_symb{iS},'color',use_color{iS},...
                        'MarkerFaceColor',use_fill{iS},'MarkerSize',use_mk_size,...
                        'linewidth',use_lw,'LineStyle','none')
                end
            end
            if options.connect_not_sig && ... % show lines connecting post hoc
                    iA~=size(use_dice,2)
                for iA2 = (iA+1):size(use_dice,2)
                    if save_FDR_p{iM}(iA,iA2) > 0.05
                        plot([iA iA2], [nanmean(use_data(:,iA)) ...
                            nanmean(use_data(:,iA2))], '-', ...
                            'color', [0.75 0.75 0.75], 'linewidth',use_lw)
                    end
                end
            end
            
            errorbar(iA,nanmean(plot_me),...
                nanstd(plot_me)/sqrt(numel(plot_me)),...
                'Marker', use_symb{iA}, 'color', use_color{iA},...
                'MarkerFaceColor', use_fill{iA}, ...
                'MarkerSize', use_mk_size,'linewidth',use_lw)
            set(gca,'XTick',use_x_tick,'XTickLabel',use_x_labels,...
                'XColor','k','YColor','k','fontsize',18)
            title(plot_q_names{iM})
            if iA == numel(analyses)
                if iM < 2
                    y_min = min(nanmean(use_dice(:,:,iM),1));
                    y_max = max(nanmean(use_dice(:,:,iM),1));
                    ylabel('Overlap (Dice coeff.)','color','k');
                    y_range = y_max - y_min;
                else
                    y_min = min(nanmean(use_minfo,1));
                    y_max = max(nanmean(use_minfo,1));
                    ylabel('Mutual info. (bits)','color','k');
                    y_range = (y_max - y_min);
                end
                                
                axis([0.5 numel(analyses)+0.5 ...
                    y_min-0.1*y_range y_max+0.1*y_range ]);
                
                patch_x = use_patch_pct_x*(numel(analyses)+0.5 - 0.5);
                patch_y = use_patch_pct_y*(y_max+0.1*y_range - y_min-0.1*y_range);
                patch([0.5-patch_x 0.5+patch_x 0.5+patch_x 0.5-patch_x],...
                    [y_min-0.1*y_range-patch_y y_min-0.1*y_range-patch_y ...
                    y_min-0.1*y_range+patch_y y_min-0.1*y_range+patch_y],...
                    'w','EdgeColor','none');
                
                set(gcf,'color','w')
                
                if iM == 2             
                    h_leg = legend(legend_names);
                    set(h_leg,'TextColor','k','FontSize',18);
                end
            end
            
        end
        
    end
    
end

%% results
results.save_FDR_p = save_FDR_p;
results.dice = dice;
results.use_dice = use_dice;
results.use_minfo = use_minfo;
results.options = options;
results.date_run = datestr(now);

end