function output = quantify_uncorr_diff( options )
% usage: output = quantify_uncorr_diff( options )
%
% mps 20210308

%% check matlab version
if datenum(version('-date')) < datenum('September 14, 2017')
    error(['This code relies on new functionality of niftiread.m, and will NOT function with versions '...
        'of Matlab older than 2017b'])
end

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

if ~isfield(options,'topDir')
    options.topDir = input('Path to data directory: ','s');
end

if ~isfield(options,'displayFigs')
    options.displayFigs = 1; % 0 = no, 1 = yes
end
if ~isfield(options,'which_error_bars')
    options.which_error_bars = 'Morey'; % Morey, sem
end

if ~isfield(options, 'connect_not_sig')
    options.connect_not_sig = 1; % gray lines connecting non-significant differences
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

output = [];

%% get subject numbers
output.subj_number = [];
for iS = 1:numel(options.subjDirs)
    output.subj_number(iS) = str2num(options.subjDirs{iS}(2:end));
end

%% load data
nSubj = numel(options.subjDirs);

analyses = {'GE_qwarp', 'GE_topup', 'fugue', 'SE_qwarp', ...
    'SE_topup'};

data_dir = fullfile(options.topDir, 'single_GE');

h_wait = waitbar(0, 'loading data, please wait...');

for iS = 1:numel(options.subjDirs)
    for iA = 1:numel(analyses)
        sub_uncorr_dir = fullfile(data_dir, options.subjDirs{iS}, options.scanSubDirs{iS}, ...
            'uncorr');
        sub_corr_dir = fullfile(data_dir, options.subjDirs{iS}, options.scanSubDirs{iS}, ...
            analyses{iA});
        
        diff_file = fullfile(sub_corr_dir, 'uncorr_min_corr.nii.gz');
        if ~exist(diff_file,'file')
            cmd = ['3dcalc -overwrite -prefix ' diff_file ' -a "' fullfile(...
                sub_uncorr_dir, ['epi_uncorr.nii.gz']) '" -b "' ...
                fullfile(sub_corr_dir, ['epi_' analyses{iA} '.nii.gz']) ...
                '" -expr "a-b"'];
            [~, result] = system(cmd);
        end
        
        load_data = round(double(niftiread(diff_file)));
        
        mask_file = fullfile(sub_corr_dir, 'uncorr_plus_corr_mask.nii.gz');
        if ~exist(mask_file,'file')
            cmd = ['3dcalc -overwrite -prefix ' mask_file ' -a "' fullfile(...
                sub_uncorr_dir, ['epi_uncorr_mask.nii.gz']) '" -b "' ...
                fullfile(sub_corr_dir, ['epi_' analyses{iA} '_mask.nii.gz']) ...
                '" -expr "step(a+b)"'];
            [~, result] = system(cmd);
        end
        
        load_mask = round(double(niftiread(mask_file)));
        mask_nan = nan(size(load_mask));
        mask_nan(load_mask == 1) = 1;
        
        load_data = load_data.*mask_nan;
        
        sub_data(iS,iA) = nanmedian( abs(load_data(:)) );
        
        waitbar(iS/nSubj, h_wait);
    end
end
close(h_wait)
%% stats
if options.displayFigs
    showFig = 'on';
else
    showFig = 'off';
end

all_data = sub_data;
all_subj = repmat([1:size(all_data,1)]',[1 size(all_data,2)]);
all_con = repmat([1:size(all_data,2)],[size(all_data,1) 1]);

[p output.anova stat] = anovan(all_data(:),{all_subj(:), all_con(:)},'random',1,...
    'model', 'full','varnames',{'subj','cond'}, 'display', showFig);

output.ttest = cell(...
    (size(all_data,2)-1), size(all_data,2) );
all_p = [];
all_iC = [];
save_FDR_p = nan((size(all_data,2)-1), size(all_data,2));
% 1 tailed t-tests for all paired comparisons
for iC1 = 1:(size(all_data,2)-1)
    for iC2 = (iC1+1):size(all_data,2)
        [h, p, ci, output.ttest{iC1,iC2}] = ...
            ttest(all_data(:, iC1), all_data(:, iC2));
        output.ttest{iC1,iC2}.p_uncorr = p;
        all_p = [all_p ; p];
        all_iC = [all_iC ; iC1 iC2];
        output.ttest{iC1,iC2}.C1_name = ...
            analyses{iC1};
        output.ttest{iC1,iC2}.C2_name = ...
            analyses{iC2};
    end
end
[sort_p, p_idx] = sort(all_p, 'ascend');
FDR_p = sort_p.*[numel(sort_p):-1:1]';
for iP = 1:numel(FDR_p)
    output.ttest{all_iC(p_idx(iP), 1), ...
        all_iC(p_idx(iP), 2)}.p_FDR = ...
        FDR_p(iP);
    save_FDR_p(all_iC(p_idx(iP), 1), ...
        all_iC(p_idx(iP), 2)) = FDR_p(iP);
end

%% plot data
use_symb = {'s','^','^','s','^'};
use_color = {'r','r',[0 0.75 0],'b','b'};
use_fill = {'r','r',[0 0.75 0],'b','b'};
use_mk_size = 10;
use_lw = 1;
use_patch_pct_y = 0.05;
use_patch_pct_x = 0.01;
use_x_tick = [1.5 3 4.5 6];
use_x_labels = {'GE oppPE','B0 FM','SE oppPE'};

legend_names = {'GE 3dQwarp', 'GE topup', 'fugue', 'SE 3dQwarp', ...
    'SE topup'};

figure
hold on

for iA = 1:numel(analyses)
    
    if strcmp(lower(options.which_error_bars),'sem') % show SEM error bars
        use_data = sub_data;
        plot_me = sub_data(:,iA);
    else strcmp(lower(options.which_error_bars),'morey') % show morey within-subj error bars
        use_data = sub_data - repmat(nanmean(sub_data,2),...
            [1 size(sub_data,2)]) + repmat(nanmean(nanmean(...
            sub_data,2),1), size(sub_data));
        plot_me = use_data(:,iA);
    end
    
    if iA == 1
        for iS = 1:numel(use_symb)
            plot(-1,-1,'Marker',use_symb{iS},'color',use_color{iS},...
                'MarkerFaceColor',use_fill{iS},'MarkerSize',use_mk_size,...
                'linewidth',use_lw,'LineStyle','none')
        end
    end
    if options.connect_not_sig && ... % show lines connecting post hoc
            iA~=size(sub_data,2)
        for iA2 = (iA+1):size(sub_data,2)
            if save_FDR_p(iA,iA2) > 0.05
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
        'MarkerSize', use_mk_size,'linewidth',use_lw,...
        'LineStyle','none')
    
end

y_min = min(nanmean(sub_data,1));
y_max = max(nanmean(sub_data,1));
y_range = y_max - y_min;
ylabel('Median abs. signal diff. (uncorr. - corr.)','color','k');
set(gca,'XTick',use_x_tick,'XTickLabel',use_x_labels,...
    'XColor','k','YColor','k','fontsize',18)

patch_x = use_patch_pct_x*(numel(analyses)+0.5 - 0.5);
patch_y = use_patch_pct_y*(y_max+0.1*y_range - y_min-0.1*y_range);
patch([0.5-patch_x 0.5+patch_x 0.5+patch_x 0.5-patch_x],...
    [y_min-0.1*y_range-patch_y y_min-0.1*y_range-patch_y ...
    y_min-0.1*y_range+patch_y y_min-0.1*y_range+patch_y],...
    'w','EdgeColor','none');

axis([0.5 numel(analyses)+0.5 ...
    y_min-0.1*y_range y_max+0.1*y_range ]);
h_leg = legend(legend_names);
set(h_leg,'TextColor','k','FontSize',18,'Location','SouthEast');
set(gcf,'color','w','POS',[360   278   700   420])


%% out
output.options = options;
end