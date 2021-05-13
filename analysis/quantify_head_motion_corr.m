function output = quantify_head_motion_corr( options )
% usage: output = quantify_head_motion_corr( options )
%
% mps 20210402

%% check matlab version
if datenum(version('-date')) < datenum('September 14, 2017')
    error(['This code relies on new functionality of niftiread.m, and will NOT function with versions '...
        'of Matlab older than 2017b'])
end

%% opt
if ~exist('options','var')
    options = [];
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

analyses = {'GE_qwarp','GE_topup','fugue','SE_qwarp','SE_topup'};

epi_files = {'epi_GE_qwarp.nii.gz', 'epi_GE_topup.nii.gz', ...
    'epi_fugue.nii.gz','epi_SE_qwarp.nii.gz','epi_SE_topup.nii.gz'};
epi_raw_files = {'PRF1_AP_gu.nii.gz', 'CSS_task3_AP_gu.nii.gz', ...
    'COP_task1_AP_gu.nii.gz'};

fm_files = {'PAGE_GE_qwarp.nii.gz', 'PAGE_GE_topup.nii.gz', ...
    'b0_mag_gu.nii.gz','APSE_SE_qwarp.nii.gz','APSE_SE_topup.nii.gz'};
fm_raw_files = {'oppPE_AP_gu.nii.gz', 'oppPE_AP_gu.nii.gz', ...
    'b0_mag_gu.nii.gz','AP_SE_gu.nii.gz','AP_SE_gu.nii.gz'};

data_dir = fullfile(options.topDir, 'separate_GE');

topup_acq_params = fullfile(data_dir,'topup_acq_params.txt');

h_wait = waitbar(0, 'loading data, please wait...');

motion_data = nan(numel(options.subjDirs), numel(analyses));

gmwm_dice_data_3epi = nan(numel(options.subjDirs), numel(analyses), 3); % 3 = 3 GE EPI scans
motion_data_3epi = nan(numel(options.subjDirs), numel(analyses), 3); % 3 = 3 GE EPI scans

for iS = 1:numel(options.subjDirs)
    for iA = 1:numel(analyses)
        subj_dir = fullfile(data_dir, options.subjDirs{iS}, options.scanSubDirs{iS});
        raw_dir = fullfile(subj_dir, 'raw');
        analysis_dir = fullfile(subj_dir, analyses{iA});
        motion_dir = fullfile(analysis_dir,'head_motion');
        
        if ~exist(motion_dir,'dir')
            mkdir(motion_dir);
        end
        
        if ~exist(fullfile(motion_dir,epi_files{iA}),'file')
            
            if ~exist(fullfile(analysis_dir,epi_files{iA}),'file')
                error(['Can''t find file: ' fullfile(analysis_dir,epi_files{iA})]);
            end
            
            cmd = ['ln -s ' fullfile(analysis_dir,epi_files{iA}) ' ' ...
                fullfile(motion_dir,epi_files{iA})];
            [~, result] = system(cmd);
        end
        
        fm_out = fullfile(motion_dir,fm_files{iA});
        fm_temp = fullfile(motion_dir,fm_raw_files{iA});
        fm_in = fullfile(raw_dir,fm_raw_files{iA});
        if ~exist(fm_out,'file')
            
            if ~exist(fm_in,'file')
                error(['Can''t find file: ' fm_in]);
            end
            
            if iA == 1 % GE qwarp
                warp_map = fullfile(analysis_dir,'blip_warp_Rev_WARP+orig'); 
                % Rev because PA data
                cmd = ['3dNwarpApply -quintic -nwarp ' warp_map ' -source '...
                    fm_in ' -prefix ' fm_out];
                [~, result] = system(cmd);
                
                cmd = ['3drefit -atrcopy ' fm_in ' IJK_TO_DICOM_REAL ' ...
                    fm_out];
                [~, result] = system(cmd);
            
            elseif iA == 2 % GE topup
                
                cmd = ['3dZeropad -overwrite -prefix ' fm_temp ' -S 1 ' ...
                    fm_in];
                [~, result] = system(cmd);
                
                topup_map = fullfile(analysis_dir,'topup_map');
                
                cmd = ['applytopup --imain=' fm_temp ' --inindex=2 --datain='...
                    topup_acq_params ' --method=jac --topup=' topup_map...
                    ' --out=' fm_out]; % inindex = 2 because PA data
                [~, result] = system(cmd);
                
                cmd = ['3dZeropad -overwrite -prefix ' fm_out ' -S -1 ' ...
                    fm_out];
                [~, result] = system(cmd);
                
            elseif iA == 3 % fugue
                
                cmd = ['ln -s ' fm_in ' ' fm_out];
                [~, result] = system(cmd);
                
            elseif iA == 4 % SE qwarp
                
                warp_map = fullfile(analysis_dir,'blip_warp_For_WARP+orig'); % For because AP data
                
                cmd = ['3dNwarpApply -quintic -nwarp ' warp_map ' -source '...
                    fm_in ' -prefix ' fm_out];
                [~, result] = system(cmd);
                
                cmd = ['3drefit -atrcopy ' fm_in ' IJK_TO_DICOM_REAL ' ...
                    fm_out];
                [~, result] = system(cmd);
                
            elseif iA == 5 % SE topup
                
                cmd = ['3dZeropad -overwrite -prefix ' fm_temp ' -S 1 ' ...
                    fm_in];
                [~, result] = system(cmd);
                
                topup_map = fullfile(analysis_dir,'topup_map');
                
                cmd = ['applytopup --imain=' fm_temp ' --inindex=1 --datain='...
                    topup_acq_params ' --method=jac --topup=' topup_map...
                    ' --out=' fm_out]; % inindex = 1 because AP data
                [~, result] = system(cmd);
                
                cmd = ['3dZeropad -overwrite -prefix ' fm_out ' -S -1 ' ...
                    fm_out];
                [~, result] = system(cmd);
                
            end
        end
        
        % do 3dvolreg to quantify motion
        motion_1D = fullfile(motion_dir,[analyses{iA} '_motion.1D']);
        if ~exist(motion_1D,'file')
            cmd = ['3dvolreg -overwrite -verbose -zpad 1 -base ' fullfile(motion_dir,...
                fm_files{iA}) ' -1Dfile ' fullfile(motion_dir,...
                ['vr_' analyses{iA} '.1D']) ' -prefix ' fullfile(motion_dir,...
                ['rm_vr_' analyses{iA} '.nii.gz']) ' -quintic '...
                '-1Dmatrix_save ' fullfile(motion_dir,['vr_mat_' ...
                analyses{iA} '.1D']) ' -maxdisp1D ' fullfile(motion_dir,...
                [analyses{iA} '_motion.1D']) ' ' fullfile(motion_dir,...
                epi_files{iA})];
            [~, result] = system(cmd);
            
        end
        cmd = ['more ' motion_1D];
        [~, result] = system(cmd);
        pattern = '# max displacement (mm) for each volume';
        idx_start = strfind(result, pattern);
        motion_data(iS,iA) = str2num(result(idx_start+length(pattern):end));
        
        %% extra analysis to check whether trend holds when applying GE or
        % SE topup correction to more distant-in-time EPI data
        if iA ~= 3 % GE or SE oppPE, not fugue
            
            % make sub folders for PRF, CSS, COP
            for iEPI = 1:numel(epi_raw_files)
                sub_motion_dir = fullfile(motion_dir,['GE_' num2str(iEPI)]);
                if ~exist(sub_motion_dir,'dir')
                    mkdir(sub_motion_dir)
                end
                
                % apply topup
                if ((iA == 1) && (iEPI == 1)) || ((iA == 2) && (iEPI == 1)) ...
                        || ((iA == 4) && (iEPI == 3)) || ((iA == 5) && (iEPI == 3))
                    % don't need to re-do dist comp, already exists, just
                    % link
                    GE_corr = fullfile(sub_motion_dir, ['epi_' analyses{iA}...
                        '_' num2str(iEPI) '.nii.gz']);
                    if ~exist(GE_corr,'file')
                        cmd = ['ln -s ' fullfile(motion_dir,['epi_' ...
                            analyses{iA} '.nii.gz']) ' ' GE_corr];
                        [~, result] = system(cmd);
                    end
                    motion_data_3epi(iS,iA,iEPI) = motion_data(iS,iA);
                    
                    % don't need to re-do segementation of GM, WM, CSF
                    gmwm = fullfile(sub_motion_dir, ['epi_' analyses{iA}...
                        '_' num2str(iEPI) '_gmwm.nii.gz']);
                    if ~exist(gmwm,'file')
                        cmd = ['ln -s ' fullfile(motion_dir,['epi_' ...
                            analyses{iA} '_gmwm.nii.gz']) ' ' gmwm];
                        [~, result] = system(cmd);
                    end
                    
                    % calculate DICE for CSF-excluded
                    dice_path = fullfile(subj_dir, 'DICE.txt');
                    if ~exist(dice_path,'file')
                        error(['can''t find file: ' dice_path]);
                    else
                        fID = fopen(dice_path,'r');
                        all_txt = textscan(fID,'%s','Delimiter','\n');
                        fclose(fID);
                        match_all = strfind(all_txt{1},[analyses{iA} ':']);
                        clear idx
                        for iM = 1:numel(match_all)
                            if match_all{iM} == 1
                                idx = iM;
                            end
                        end
                        gmwm_dice_data_3epi(iS,iA,iEPI) = str2num(all_txt{1}{idx+2}); % gmwm dice
                    end
                    
                else
                    % do dist comp using new EPI
                    fm_out = fullfile(sub_motion_dir,['epi_' analyses{iA}...
                        '_' num2str(iEPI) '.nii.gz']);
                    fm_temp = fullfile(sub_motion_dir,epi_raw_files{iEPI});
                    fm_in = fullfile(raw_dir,epi_raw_files{iEPI});
                    
                    if ~exist(fm_out, 'file')
                        if iA == 1 || iA == 4 % GE or SE 3dQwarp
                            warp_map = fullfile(analysis_dir,'blip_warp_For_WARP+orig');
                            % For because AP data
                            cmd = ['3dNwarpApply -quintic -nwarp ' warp_map ' -source '...
                                fm_in ' -prefix ' fm_out];
                            [~, result] = system(cmd);
                            
                            cmd = ['3drefit -atrcopy ' fm_in ' IJK_TO_DICOM_REAL ' ...
                                fm_out];
                            [~, result] = system(cmd);
                            
                        elseif iA == 2 || iA == 5 % GE or SE topup
                            cmd = ['3dZeropad -overwrite -prefix ' fm_temp ' -S 1 ' ...
                                fm_in];
                            [~, result] = system(cmd);

                            topup_map = fullfile(analysis_dir,'topup_map');

                            cmd = ['applytopup --imain=' fm_temp ' --inindex=1 --datain='...
                                topup_acq_params ' --method=jac --topup=' topup_map...
                                ' --out=' fm_out]; % inindex = 1 because AP data
                            [~, result] = system(cmd);

                            cmd = ['3dZeropad -overwrite -prefix ' fm_out ' -S -1 ' ...
                                fm_out];
                            [~, result] = system(cmd);
                        end
                    end
                
                    % align T1
                    if ~exist(fullfile(sub_motion_dir, ['epi_' analyses{iA}...
                            '_' num2str(iEPI) '_al_3T_mat.aff12.1D']),'file')
                        align_script = fullfile(options.gitDir,'align_one.sh');
                        out_txt = fullfile(sub_motion_dir, 'out.txt');
                        if exist(out_txt,'file')
                            delete(out_txt);
                        end
                        
                        cmd = ['touch ' out_txt];
                        [~, result] = system(cmd);
                        
                        cmd = ['chmod a+w ' out_txt];
                        [~, result] = system(cmd);
                        
                        cmd = ['tcsh ' align_script ' ' subj_dir ' ' analyses{iA} ...
                            ' ' num2str(iEPI) ' >> ' out_txt];
                        [~, result] = system(cmd);
                    end
                    
                    % do segementation of GM, WM, CSF -- adapt segment_all_six.sh
                    segment_script = fullfile(options.gitDir,'segment_one.sh');

                    if ~exist(fullfile(sub_motion_dir, ['epi_' analyses{iA}...
                            '_' num2str(iEPI) '_seg.nii.gz']),'file')
                        cmd = ['tcsh ' segment_script ' ' subj_dir ' ' analyses{iA} ...
                            ' ' num2str(iEPI) ' >> ' out_txt];
                        [~, result] = system(cmd);
                    end
                    
                    % calculate DICE for CSF-excluded
                    new_dice_file = fullfile(sub_motion_dir, 'DICE.txt');
                    if ~exist(new_dice_file,'file')
                        error(['can''t find file: ' new_dice_file]);
                    else
                        fID = fopen(new_dice_file,'r');
                        all_txt = textscan(fID,'%s','Delimiter','\n');
                        fclose(fID);
                        match_all = strfind(all_txt{1},[analyses{iA} ':']);
                        clear idx
                        for iM = 1:numel(match_all)
                            if match_all{iM} == 1
                                idx = iM;
                            end
                        end
                        gmwm_dice_data_3epi(iS,iA,iEPI) = str2num(all_txt{1}{idx+2}); % gmwm dice
                    end
                    
                    % do 3dvolreg on all 3 GE AP scans, with corrected PA GE or AP SE as the base
                    motion_1D = fullfile(sub_motion_dir,[analyses{iA} '_' ...
                        num2str(iEPI) '_motion.1D']);
                    if ~exist(motion_1D,'file')
                        fm_out_vr = fullfile(sub_motion_dir,['epi_' analyses{iA}...
                            '_' num2str(iEPI) '_vr.nii.gz']);
                        cmd = ['3dvolreg -overwrite -verbose -zpad 1 -base ' ...
                            fullfile(motion_dir,fm_files{iA}) ' -1Dfile ' ...
                            fullfile(sub_motion_dir,...
                            ['vr_' analyses{iA} '_' num2str(iEPI) '.1D']) ...
                            ' -prefix ' fm_out_vr ' -quintic -1Dmatrix_save ' ...
                            fullfile(sub_motion_dir,['vr_mat_' analyses{iA} ...
                            '_' num2str(iEPI) '.1D']) ' -maxdisp1D ' ...
                            motion_1D ' ' fm_out];
                        [~, result] = system(cmd);
                        
                    end
                    cmd = ['more ' motion_1D];
                    [~, result] = system(cmd);
                    pattern = '# max displacement (mm) for each volume';
                    idx_start = strfind(result, pattern);
                    motion_data_3epi(iS,iA,iEPI) = str2num(result(idx_start+length(pattern):end));
                end
            end
        end
    end
    
waitbar(iS/nSubj, h_wait);
end
close(h_wait)

%% run main distortion comp analysis, pull in Dice and Mutual Info
dc_opt = [];
dc_opt.displayFigs = 0;
dc_opt.subjDirs = options.subjDirs;
dc_opt.scanSubDirs = options.scanSubDirs;
dc_opt.topDir = options.topDir;
dc_opt.gitDir = options.gitDir;

dc_data = dist_comp_analysis( dc_opt );

%% stats
if options.displayFigs
    showFig = 'on';
else
    showFig = 'off';
end

all_data = motion_data;
all_subj = repmat([1:size(all_data,1)]',[1 size(all_data,2)]);
all_con = repmat([1:size(all_data,2)],[size(all_data,1) 1]);

[p output.anova_motion stat] = anovan(all_data(:),{all_subj(:), all_con(:)},'random',1,...
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

% stats for 3 EPI x 2 (GE & SE oppPE) x 2 (3dQwarp & topup) - motion
all_data = motion_data_3epi(:,[1 2 4 5],:); % 1 = GE 3dQwarp, 2 = GE topup, 4 = SE 3dQwarp, 5 = SE topup
all_subj = repmat([1:size(all_data,1)]',[1 size(all_data,2) size(all_data,3)]);
all_GE_SE = repmat([1 1 2 2],[size(all_data,1) 1  size(all_data,3)]);
all_qwarp_topup = repmat([1 2 1 2],[size(all_data,1) 1  size(all_data,3)]);
all_scan = cat(3, ones(size(all_data,1),size(all_data,2)), ...
    2*ones(size(all_data,1),size(all_data,2)), ...
    3*ones(size(all_data,1),size(all_data,2)));

[p output.anova_motion_3epi stat] = anovan(all_data(:),{all_subj(:), ...
    all_GE_SE(:), all_qwarp_topup(:), all_scan(:)},'random', 1,...
    'model', 'full','varnames',{'subj','GE SE', '3dQwarp topup', 'epi scan'},...
    'display', showFig);


% stats for 3 EPI x 2 (GE & SE oppPE) - Dice
all_data = gmwm_dice_data_3epi(:,[1 2 4 5],:); % 2 = GE topup, 5 = SE topup
all_subj = repmat([1:size(all_data,1)]',[1 size(all_data,2) size(all_data,3)]);
all_GE_SE = repmat([1 1 2 2],[size(all_data,1) 1  size(all_data,3)]);
all_qwarp_topup = repmat([1 2 1 2],[size(all_data,1) 1  size(all_data,3)]);
all_scan = cat(3, ones(size(all_data,1),size(all_data,2)), ...
    2*ones(size(all_data,1),size(all_data,2)), ...
    3*ones(size(all_data,1),size(all_data,2)));

[p output.anova_dice_3epi stat] = anovan(all_data(:),{all_subj(:), ...
    all_GE_SE(:), all_qwarp_topup(:), all_scan(:)},'random', 1,...
    'model', 'full','varnames',{'subj','GE SE', '3dQwarp topup', 'epi scan'},...
    'display', showFig);

%% compute correlations and plot data

subplot_vals2 = {[0.09 0.6 0.39 0.34],...
    [0.59 0.6 0.39 0.34],...
    [0.09 0.1 0.39 0.34]};

use_mk_size = 10;
use_lw = 1;
use_patch_pct_y = 0.05;
use_patch_pct_x = 0.01;

use_symb = {'s','^','^','s','^'};
use_color = {'r','r',[0 0.75 0],'b','b'};
use_fill = {[1 0.5 0.5],[1 0.5 0.5],[0.5 0.75 0.5],[0.5 0.5 1],[0.5 0.5 1]};

plot_q_names = {'Dice (whole-brain masks) vs. motion','Dice (CSF-excl. masks) vs. motion',...
    'Mutual info. vs. motion'};

dc_names = {'GE_qwarp','GE_topup','fugue','SE_qwarp','SE_topup'};
legend_names = {'GE 3dQwarp', 'GE topup', 'fugue','SE 3dQwarp', ...
    'SE topup'};

cond_idx = [1 2 3]; % mask, GMWM, mututal info

% compute correlation between head motion and distortion comp metrics
for iC = 1:numel(cond_idx)
    all_x = motion_data;
    if iC < 3
        all_y = squeeze(dc_data.use_dice(:,1:5,iC));
    else
        all_y = dc_data.use_minfo(:,1:5);
    end
    [output.corr.motion.r(iC), output.corr.motion.p(iC)] = corr(all_x(:), all_y(:),...
        'type','spearman');
    output.corr.motion.df(iC) = numel(all_x)-2;
    
    [poly_fit(iC,:)] = polyfit(all_x(:), all_y(:), 1);
    
end

% compute correlation between head motion and distortion comp metrics for 3
% epi data
all_x = motion_data_3epi(:,[1 2 4 5],:);
all_y = gmwm_dice_data_3epi(:,[1 2 4 5],:);

[output.corr.epi3.r, output.corr.epi3.p] = corr(all_x(:), all_y(:),...
    'type','spearman');
output.corr.epi3.df = numel(all_x)-2;

[poly_fit_3epi] = polyfit(all_x(:), all_y(:), 1);

y_min = [];
y_max = [];

if options.displayFigs
    h2 = figure;
    
    for iA = 1:numel(dc_names)
        for iM = cond_idx
            
            if iM < 3
                plot_y = squeeze(dc_data.use_dice(:,iA,iM));
                all_y = dc_data.use_dice(:,:,iM);
                y_min = min(all_y(:));
                y_max = max(all_y(:));
            else
                plot_y = dc_data.use_minfo(:,iA);
                all_y = dc_data.use_minfo;
                y_min = min(all_y(:));
                y_max = max(all_y(:));
            end
            
            plot_x = motion_data(:,iA);
            
            subplot('position',subplot_vals2{iM}); hold on
            if iM == 3 && iA == 1
                for iS = 1:numel(use_symb)
                    plot(-100,-100,'Marker',use_symb{iS},'color',use_color{iS},...
                        'MarkerFaceColor',use_fill{iS},'MarkerSize',use_mk_size,...
                        'linewidth',use_lw,'LineStyle','none')
                end
            end
            
            plot(plot_x,plot_y,...
                'Marker', use_symb{iA}, 'color', use_color{iA},...
                'MarkerFaceColor', use_fill{iA}, ...
                'MarkerSize', use_mk_size,'linewidth',use_lw,...
                'LineStyle','none')
            set(gca,'XColor','k','YColor','k','fontsize',18)
            title(plot_q_names{iM})
            if iA == numel(dc_names)
                
                if iM < 3
                    ylabel('Overlap (Dice coeff.)','color','k');
                else
                    ylabel('Mutual info. (bits)','color','k');
                end
                xlabel('Residual head motion (mm)','color','k')
                
                y_range = y_max - y_min;
                
                x_min = min(motion_data(:));
                x_max = max(motion_data(:));
                x_range = x_max - x_min;
                
                fit_x = [x_min x_max];
                fit_y = poly_fit(iM,1).*fit_x + poly_fit(iM,2);
                
                plot(fit_x, fit_y, 'k-', 'linewidth', 2)
                
                text(fit_x(2) - x_range*0.25,y_min + y_range*0.125,...
                    ['r = ' num2str(round(output.corr.motion.r(iM)*100)/100)],'fontsize',18)
                text(fit_x(2) - x_range*0.25,y_min,...
                    ['p = ' num2str(round(output.corr.motion.p(iM)*1000)/1000)],'fontsize',18)
                
                axis([x_min-0.1*x_range x_max+0.1*x_range ...
                    y_min-0.1*y_range y_max+0.1*y_range ]);
                
                patch_x = use_patch_pct_x*(x_max+0.1*x_range - x_min-0.1*x_range);
                patch_y = use_patch_pct_y*(y_max+0.1*y_range - y_min-0.1*y_range);
                patch([x_min-0.1*x_range-patch_x x_min-0.1*x_range+patch_x ...
                    x_min-0.1*x_range+patch_x x_min-0.1*x_range-patch_x],...
                    [y_min-0.1*y_range-patch_y y_min-0.1*y_range-patch_y ...
                    y_min-0.1*y_range+patch_y y_min-0.1*y_range+patch_y],...
                    'w','EdgeColor','none');
                
                set(gcf,'color','w')
                
                if iM == 3
                    h_leg = legend(legend_names);
                    set(h_leg,'TextColor','k','FontSize',18);
                end
            end
            
        end
        
    end
    
    %% plot 3 epi motion
    figure;
    use_analyses = [1 2 4 5];
    for iPlot = 1:numel(use_analyses)
        subplot(2,2,iPlot);
        hold on
        
        for iEPI = 1:size(motion_data_3epi,3)
            if strcmp(lower(options.which_error_bars),'sem') % show SEM error bars
                use_data = motion_data_3epi(:,use_analyses(iPlot),:);
                plot_me = motion_data_3epi(:,use_analyses(iPlot),iEPI);
            else strcmp(lower(options.which_error_bars),'morey') % show morey within-subj error bars
                use_data = motion_data_3epi - ...
                    repmat(nanmean(nanmean(motion_data_3epi,3),2),...
                    [1 size(motion_data_3epi,2) size(motion_data_3epi,3)]) + ...
                    repmat(nanmean(nanmean(nanmean(motion_data_3epi,3),2),1), ...
                    [size(motion_data_3epi,1) size(motion_data_3epi,2) ...
                    size(motion_data_3epi,3)]); % subtract individual subject mean, add group mean
                plot_me = use_data(:,use_analyses(iPlot),iEPI);
            end
            errorbar(iEPI,nanmean(plot_me),...
                nanstd(plot_me)/sqrt(numel(plot_me)),...
                'Marker', use_symb{use_analyses(iPlot)}, ...
                'color', use_color{use_analyses(iPlot)},...
                'MarkerFaceColor', use_fill{use_analyses(iPlot)}, ...
                'MarkerSize', use_mk_size,'linewidth',use_lw)
        end
        
        box off
        set(gca,'XTick',[1:3],'XColor','k','YColor','k','fontsize',18)
        xlabel('EPI scan','color','k')
        ylabel('Residual head motion','color','k')
        title(legend_names{use_analyses(iPlot)},'color','k')
        axis([0.5 3.5 -0.25 5]) 
    end
    set(gcf,'color','w')
        
    %% plot 3 epi Dice
    figure;
    use_analyses = [1 2 4 5];
    for iPlot = 1:numel(use_analyses)
        subplot(2,2,iPlot);
        hold on
        
        for iEPI = 1:size(gmwm_dice_data_3epi,3)
            if strcmp(lower(options.which_error_bars),'sem') % show SEM error bars
                use_data = gmwm_dice_data_3epi(:,use_analyses(iPlot),:);
                plot_me = gmwm_dice_data_3epi(:,use_analyses(iPlot),iEPI);
            else strcmp(lower(options.which_error_bars),'morey') % show morey within-subj error bars
                use_data = gmwm_dice_data_3epi - ...
                    repmat(nanmean(nanmean(gmwm_dice_data_3epi,3),2),...
                    [1 size(gmwm_dice_data_3epi,2) size(gmwm_dice_data_3epi,3)]) + ...
                    repmat(nanmean(nanmean(nanmean(gmwm_dice_data_3epi,3),2),1), ...
                    [size(gmwm_dice_data_3epi,1) size(gmwm_dice_data_3epi,2) ...
                    size(gmwm_dice_data_3epi,3)]); % subtract individual subject mean, add group mean
                plot_me = use_data(:,use_analyses(iPlot),iEPI);
            end
            errorbar(iEPI,nanmean(plot_me),...
                nanstd(plot_me)/sqrt(numel(plot_me)),...
                'Marker', use_symb{use_analyses(iPlot)}, ...
                'color', use_color{use_analyses(iPlot)},...
                'MarkerFaceColor', use_fill{use_analyses(iPlot)}, ...
                'MarkerSize', use_mk_size,'linewidth',use_lw)
        end
        
        box off
        set(gca,'XTick',[1:3],'XColor','k','YColor','k','fontsize',18)
        xlabel('EPI scan','color','k')
        ylabel('Dice (CSF-excl. masks)','color','k')
        title(legend_names{use_analyses(iPlot)},'color','k')
        axis([0.5 3.5 0.865 0.885]) 
    end
    set(gcf,'color','w')
    
    %% plot correlation of 3 epi motion & dice
    figure; 
    subplot(2,1,1); hold on
    for iA = [1 2 4 5]
        plot_x = motion_data_3epi(:,iA,:);
        plot_y = gmwm_dice_data_3epi(:,iA,:);
        
        plot(plot_x(:), plot_y(:), ...
            'Marker', use_symb{iA}, 'color', use_color{iA},...
            'MarkerFaceColor', use_fill{iA}, ...
            'MarkerSize', use_mk_size,'linewidth',use_lw,...
            'LineStyle','none')
    end
    y_min = min(gmwm_dice_data_3epi(:));
    y_max = max(gmwm_dice_data_3epi(:));
    y_range = y_max - y_min;
    
    x_min = min(motion_data_3epi(:));
    x_max = max(motion_data_3epi(:));
    x_range = x_max - x_min;
    
    fit_x = [x_min x_max];
    fit_y = poly_fit_3epi(1).*fit_x + poly_fit_3epi(2);
    
    plot(fit_x, fit_y, 'k-', 'linewidth', 2)
    
    text(fit_x(2) - x_range*0.35,y_max + y_range*0.125,...
        ['r = ' num2str(round(output.corr.epi3.r*100)/100)],'fontsize',18)
    text(fit_x(2) - x_range*0.35,y_max,...
        ['p = ' num2str(round(output.corr.epi3.p*10000)/10000)],'fontsize',18)
    
    set(gca,'XColor','k','YColor','k','fontsize',18)
    
    axis([x_min-0.1*x_range x_max+0.1*x_range ...
                    y_min-0.1*y_range y_max+0.1*y_range ]);
    
    patch_x = use_patch_pct_x*(x_max+0.1*x_range - x_min-0.1*x_range);
    patch_y = use_patch_pct_y*(y_max+0.1*y_range - y_min-0.1*y_range);
    patch([x_min-0.1*x_range-patch_x x_min-0.1*x_range+patch_x ...
        x_min-0.1*x_range+patch_x x_min-0.1*x_range-patch_x],...
        [y_min-0.1*y_range-patch_y y_min-0.1*y_range-patch_y ...
        y_min-0.1*y_range+patch_y y_min-0.1*y_range+patch_y],...
        'w','EdgeColor','none');
    
    set(gcf,'color','w')
    
     ylabel('Dice (CSF-excl. masks)','color','k');
     xlabel('Residual head motion (mm)','color','k')
     
     %% plot mean dice & motion for GE & SE across 3 epi, plus predicted
     
     subplot(2,1,2); hold on
     
     use_color2 = {[1 0.85 0.85],[1 0.85 0.85],[0.5 0.75 0.5],[0.85 0.85 1],[0.85 0.85 1]};
     use_fill2 = {[1 0.85 0.85],[1 0.85 0.85],[0.5 0.75 0.5],[0.85 0.85 1],[0.85 0.85 1]};
     
     % re-plot all points
     for iA = [1 2 4 5]
         plot_x = motion_data_3epi(:,iA,:);
         plot_y = gmwm_dice_data_3epi(:,iA,:);
         
         plot(plot_x(:), plot_y(:), ...
             'Marker', use_symb{iA}, 'color', use_color2{iA},...
             'MarkerFaceColor', use_fill2{iA}, ...
             'MarkerSize', use_mk_size,'linewidth',use_lw,...
             'LineStyle','none')
     end
     
     colors = {[0.75 0 0],[0 0 0.75]};
     colors2 = {'r','b'};
     
     % correlate motion & dice for GE and SE separately
     GE_SE_idx = {[1 2],[4 5]};
     for iA = 1:numel(GE_SE_idx)
         
         all_x = motion_data_3epi(:,GE_SE_idx{iA},:);
         all_y = gmwm_dice_data_3epi(:,GE_SE_idx{iA},:);
         
         [output.corr.epi3_GE_SE.r(iA), output.corr.epi3_GE_SE.p(iA)] = corr(all_x(:), all_y(:),...
             'type','spearman');
         output.corr.epi3_GE_SE.df(iA) = numel(all_x)-2;
         
         [poly_fit_3epi_GE_SE(iA,:)] = polyfit(all_x(:), all_y(:), 1);
         
         y_min_GE_SE = min(all_y(:));
         y_max_GE_SE = max(all_y(:));
         y_range_GE_SE = y_max_GE_SE - y_min_GE_SE;
         
         x_min_GE_SE = min(all_x(:));
         x_max_GE_SE = max(all_x(:));
         x_range_GE_SE = x_max_GE_SE - x_min_GE_SE;
         
         fit_x = [x_min_GE_SE x_max_GE_SE];
         fit_y = poly_fit_3epi_GE_SE(iA,1).*fit_x + poly_fit_3epi_GE_SE(iA,2);
         
         plot(fit_x, fit_y, '-', 'color', colors{iA}, 'linewidth', 2)
         
         % now plot means & error bars for x (motion) and y (dice)
         if strcmp(lower(options.which_error_bars),'sem') % show SEM error bars
             use_data = motion_data_3epi;
             plot_x = use_data(:,GE_SE_idx{iA},:);
             
             use_data = gmwm_dice_data_3epi;
             plot_y = use_data(:,GE_SE_idx{iA},:);
         else strcmp(lower(options.which_error_bars),'morey') % show morey within-subj error bars
             use_data = motion_data_3epi - ...
                 repmat(nanmean(nanmean(motion_data_3epi,3),2),...
                 [1 size(motion_data_3epi,2) size(motion_data_3epi,3)]) + ...
                 repmat(nanmean(nanmean(nanmean(motion_data_3epi,3),2),1), ...
                 [size(motion_data_3epi,1) size(motion_data_3epi,2) ...
                 size(motion_data_3epi,3)]); % subtract individual subject mean, add group mean
             plot_x = use_data(:,GE_SE_idx{iA},:);
             
             use_data = gmwm_dice_data_3epi - ...
                 repmat(nanmean(nanmean(gmwm_dice_data_3epi,3),2),...
                 [1 size(gmwm_dice_data_3epi,2) size(gmwm_dice_data_3epi,3)]) + ...
                 repmat(nanmean(nanmean(nanmean(gmwm_dice_data_3epi,3),2),1), ...
                 [size(gmwm_dice_data_3epi,1) size(gmwm_dice_data_3epi,2) ...
                 size(gmwm_dice_data_3epi,3)]); % subtract individual subject mean, add group mean
             plot_y = use_data(:,GE_SE_idx{iA},:);
             
         end
         
         x_eb = [nanmean(plot_x(:)) - nanstd(plot_x(:))/sqrt(size(motion_data_3epi,1)) ...
             nanmean(plot_x(:)) + nanstd(plot_x(:))/sqrt(size(motion_data_3epi,1))];
         y_eb = [nanmean(plot_y(:)) - nanstd(plot_y(:))/sqrt(size(motion_data_3epi,1)) ...
             nanmean(plot_y(:)) + nanstd(plot_y(:))/sqrt(size(motion_data_3epi,1))];
         
         plot(x_eb, repmat(nanmean(plot_y(:)),[1 2]), '-', 'color', colors2{iA}, ...
             'linewidth', 2)
         plot(repmat(nanmean(plot_x(:)),[1 2]), y_eb, '-', 'color', colors2{iA}, ...
             'linewidth', 2)
         plot(nanmean(plot_x(:)), nanmean(plot_y(:)), 'o', 'color', colors2{iA}, ...
             'MarkerFaceColor', colors2{iA})
         
     end
    set(gca,'XColor','k','YColor','k','fontsize',18)
    
    y_min = 0.865 + 0.1*y_range; % hard-code this to get a zoomed-in y-axis
    y_max = 0.885 - 0.1*y_range;
    
    y_range = y_max - y_min;
    
    y_min = 0.865 + 0.1*y_range; % hard-code this to get a zoomed-in y-axis
    y_max = 0.885 - 0.1*y_range;
    
    y_range = y_max - y_min;
    
    axis([x_min-0.1*x_range x_max+0.1*x_range ...
                    y_min-0.1*y_range y_max+0.1*y_range ]);
                
    y_range = y_max - y_min;
                    
    patch_x = use_patch_pct_x*(x_max+0.1*x_range - x_min-0.1*x_range);
    patch_y = use_patch_pct_y*(y_max+0.1*y_range - y_min-0.1*y_range);
    patch([x_min-0.1*x_range-patch_x x_min-0.1*x_range+patch_x ...
        x_min-0.1*x_range+patch_x x_min-0.1*x_range-patch_x],...
        [y_min-0.1*y_range-patch_y y_min-0.1*y_range-patch_y ...
        y_min-0.1*y_range+patch_y y_min-0.1*y_range+patch_y],...
        'w','EdgeColor','none');
    
    set(gcf,'color','w','POS',[360    30   465   665])
    
    ylabel('Dice (CSF-excl. masks)','color','k');
    xlabel('Residual head motion (mm)','color','k')
     
end

%% out
output.options = options;
end