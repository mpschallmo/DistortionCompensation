#!/bin/csh
if ( $#argv > 0 ) then
    set subj = $argv[1]
else
    echo "Please specify subj as an input variable"
endif

set topDir = '/*****PATH TO ANALYSIS DIRECTORY*****/HCP/' # n.b. this must be set before running!
set subjDir = {$topDir}{$subj}'/'
set procDir = {$subjDir}'proc/'
set rawDir = {$procDir}'raw/'

set topup_acq_params = '/*****PATH TO GIT REPO*****/ref_files/topup_acq_params.txt' # n.b. this must be set before running!

mkdir {$procDir}
mkdir {$rawDir}
cd {$rawDir}

# set directories
set epi_uncorr_dir = {$subjDir}{$subj}_7T_tfMRI_RETBAR1_unproc/{$subj}/unprocessed/7T/tfMRI_RETBAR1_AP/
set T1_uncorr_dir = {$subjDir}structural_extended/{$subj}/T1w/{$subj}/mri/orig/

# get EPI
# pull out first 3 TRs
set epi_all = {$epi_uncorr_dir}{$subj}_7T_tfMRI_RETBAR1_AP.nii.gz
set epi = {$rawDir}epi_uncorr.nii.gz

# take the temporal median, put data into raw dir
3dTstat -overwrite -median -prefix {$epi} {$epi_all}'[0..2]'

# get SE
set SE_AP_all = {$epi_uncorr_dir}{$subj}_7T_SpinEchoFieldMap_AP.nii.gz
set SE_PA_all = {$epi_uncorr_dir}{$subj}_7T_SpinEchoFieldMap_PA.nii.gz
set SE_AP = {$rawDir}SE_AP.nii.gz
set SE_PA = {$rawDir}SE_PA.nii.gz

# take the temporal median, put data into raw dir
3dTstat -overwrite -median -prefix {$SE_AP} {$SE_AP_all}'[0..2]'
3dTstat -overwrite -median -prefix {$SE_PA} {$SE_PA_all}'[0..2]'

# get T1
set T1 = {$T1_uncorr_dir}{$subj}_T1_orig.nii.gz

# remove skull, move orig T1 to raw dir
set T1_strip = {$rawDir}T1_orig_strip.nii.gz
3dSkullStrip -orig_vol -input {$T1} -prefix {$T1_strip}

# remove spatial inhomogeneity in T1
set T1_uni = {$rawDir}T1_orig_strip_uni.nii.gz
3dUnifize -overwrite -prefix {$T1_uni} {$T1_strip}

# #####
# first do 6 param uncorr
set analysis = uncorr
set useDir = {$procDir}{$analysis}'/'
mkdir $useDir
cd $useDir

# copy files
set use_epi = {$useDir}epi_{$analysis}.nii.gz
set use_T1 = {$useDir}T1_strip_uni.nii.gz
3dcopy {$epi} {$use_epi}
3dcopy {$T1_uni} {$use_T1}

# align unproc T1 and epi
align_epi_anat.py -epi2anat           \
     -rigid_body                      \
     -anat {$use_T1}                  \
     -suffix _al_3T                   \
     -epi {$use_epi}                  \
     -epi_base 0                      \
     -epi_strip 3dAutomask            \
     -ginormous_move                  \
     -volreg off                      \
     -tshift off                      \
     -anat_has_skull no               \
     -overwrite
# do 6 param, to match our other analyses

set xform = {$useDir}epi_{$analysis}_al_3T_mat.aff12.1D
set ixform = {$useDir}inv_epi_{$analysis}_al_3T_mat.aff12.1D

cat_matvec -ONELINE {$xform} -I > $ixform

set T1_al = {$useDir}T1_{$analysis}.nii.gz

# resample 3T anat to match 7T EPI
3dAllineate                   \
    -overwrite                \
    -source {$use_T1}         \
    -base {$use_epi}          \
    -1Dmatrix_apply {$ixform} \
    -prefix {$T1_al}          \

# mask uncorr T1
set T1_mask = {$useDir}T1_{$analysis}_mask.nii.gz
3dAutomask -overwrite  \
    -prefix {$T1_mask} \
    {$T1_al}

# mask uncorr epi
set epi_mask = {$useDir}epi_{$analysis}_mask.nii.gz
3dAutomask -overwrite -prefix {$epi_mask} \
    {$use_epi}


####
# then do 12 param uncorr
set analysis = 12param
set useDir = {$procDir}{$analysis}'/'
mkdir $useDir
cd $useDir

# copy files
set use_epi = {$useDir}epi_{$analysis}.nii.gz
set use_T1 = {$useDir}T1_strip_uni.nii.gz
3dcopy {$epi} {$use_epi}
3dcopy {$T1_uni} {$use_T1}

# align unproc T1 and epi
align_epi_anat.py -epi2anat           \
     -anat {$use_T1}                  \
     -suffix _al_3T                   \
     -epi {$use_epi}                  \
     -epi_base 0                      \
     -epi_strip 3dAutomask            \
     -ginormous_move                  \
     -volreg off                      \
     -tshift off                      \
     -anat_has_skull no               \
     -overwrite
# do 12 param, to match what HCP did

set xform = {$useDir}epi_{$analysis}_al_3T_mat.aff12.1D
set ixform = {$useDir}inv_epi_{$analysis}_al_3T_mat.aff12.1D

cat_matvec -ONELINE {$xform} -I > $ixform

set T1_al = {$useDir}T1_{$analysis}.nii.gz

# resample 3T anat to match 7T EPI
3dAllineate                   \
    -overwrite                \
    -source {$use_T1}         \
    -base {$use_epi}          \
    -1Dmatrix_apply {$ixform} \
    -prefix {$T1_al}          \

# mask 12param T1
set T1_mask = {$useDir}T1_{$analysis}_mask.nii.gz
3dAutomask -overwrite  \
    -prefix {$T1_mask} \
    {$T1_al}

# mask 12param epi
set epi_mask = {$useDir}epi_{$analysis}_mask.nii.gz
3dAutomask -overwrite -prefix {$epi_mask} \
    {$use_epi}



####
# do SE_topup
set analysis = SE_topup
set useDir = {$procDir}{$analysis}'/'
mkdir $useDir
cd $useDir

# topup needs even number of slices for subsampling
set use_epi = {$useDir}epi_{$analysis}.nii.gz
set use_T1 = {$useDir}T1_strip_uni.nii.gz
3dcopy {$T1_uni} {$use_T1}

set fwd_med = {$SE_AP}
set rev_med = {$SE_PA}
set applyTo = {$epi}
set padded_epi = {$useDir}padded_epi.nii.gz
3dZeropad -overwrite -prefix $fwd_med -S 1 $fwd_med
3dZeropad -overwrite -prefix $rev_med -S 1 $rev_med
3dZeropad -overwrite -prefix $padded_epi -S 1 $applyTo

set both_vol = {$useDir}both_volumes.nii.gz
3dTcat -overwrite -prefix $both_vol $fwd_med $rev_med

set topup_map = {$useDir}topup_map
topup --imain=$both_vol --datain=$topup_acq_params --config=b02b0.cnf \
    --out=$topup_map
applytopup --imain=$padded_epi --inindex=1 --datain=$topup_acq_params \
    --method=jac --topup=$topup_map --out={$useDir}padded_topped_up
3dZeropad -overwrite -prefix {$use_epi} -S -1 \
    {$useDir}padded_topped_up.nii.gz

# get rid of any negative numbers (topup makes them)
3dcalc -overwrite -prefix {$use_epi} -a {$use_epi} -expr 'step(a)*a'

# align unproc T1 and epi
align_epi_anat.py -epi2anat           \
     -rigid_body                      \
     -anat {$use_T1}                  \
     -suffix _al_3T                   \
     -epi {$use_epi}                  \
     -epi_base 0                      \
     -epi_strip 3dAutomask            \
     -ginormous_move                  \
     -volreg off                      \
     -tshift off                      \
     -anat_has_skull no               \
     -overwrite
# do 6 param, to match our other analyses

set xform = {$useDir}epi_{$analysis}_al_3T_mat.aff12.1D
set ixform = {$useDir}inv_epi_{$analysis}_al_3T_mat.aff12.1D

cat_matvec -ONELINE {$xform} -I > $ixform

set T1_al = {$useDir}T1_{$analysis}.nii.gz

# resample 3T anat to match 7T EPI
3dAllineate                   \
    -overwrite                \
    -source {$use_T1}         \
    -base {$use_epi}          \
    -1Dmatrix_apply {$ixform} \
    -prefix {$T1_al}          \

# mask uncorr T1
set T1_mask = {$useDir}T1_{$analysis}_mask.nii.gz
3dAutomask -overwrite  \
    -prefix {$T1_mask} \
    {$T1_al}

# mask uncorr epi
set epi_mask = {$useDir}epi_{$analysis}_mask.nii.gz
3dAutomask -overwrite -prefix {$epi_mask} \
    {$use_epi}


####
# do SE qwarp
set analysis = SE_qwarp
set useDir = {$procDir}{$analysis}'/'
mkdir $useDir
cd $useDir

set use_epi = {$useDir}epi_{$analysis}.nii.gz
set use_T1 = {$useDir}T1_strip_uni.nii.gz
3dcopy {$T1_uni} {$use_T1}

set fwd_med = {$SE_AP}
set rev_med = {$SE_PA}
set applyTo = {$epi}

# automask the median datasets 
3dAutomask -apply_prefix {$useDir}/fwd_masked.nii.gz $fwd_med
3dAutomask -apply_prefix {$useDir}/rev_masked.nii.gz $rev_med

# compute the midpoint warp between the median datasets
3dQwarp -plusminus -pmNAMES Rev For                       \
    -pblur 0.05 0.05 -blur -1 -1                          \
    -noweight -minpatch 9                                 \
    -source {$useDir}/rev_masked.nii.gz                   \
    -base   {$useDir}/fwd_masked.nii.gz                   \
    -prefix {$useDir}/blip_warp

3dNwarpApply -quintic -nwarp {$useDir}/blip_warp_For_WARP+orig      \
         -source {$epi}                                             \
         -prefix {$use_epi}
         # mps 20191009
3drefit -atrcopy {$epi} IJK_TO_DICOM_REAL      \
             {$use_epi}
             # mps 20191009

# align unproc T1 and epi
align_epi_anat.py -epi2anat           \
     -rigid_body                      \
     -anat {$use_T1}                  \
     -suffix _al_3T                   \
     -epi {$use_epi}                  \
     -epi_base 0                      \
     -epi_strip 3dAutomask            \
     -ginormous_move                  \
     -volreg off                      \
     -tshift off                      \
     -anat_has_skull no               \
     -overwrite
# do 6 param, to match our other analyses

set xform = {$useDir}epi_{$analysis}_al_3T_mat.aff12.1D
set ixform = {$useDir}inv_epi_{$analysis}_al_3T_mat.aff12.1D

cat_matvec -ONELINE {$xform} -I > $ixform

set T1_al = {$useDir}T1_{$analysis}.nii.gz

# resample 3T anat to match 7T EPI
3dAllineate                   \
    -overwrite                \
    -source {$use_T1}         \
    -base {$use_epi}          \
    -1Dmatrix_apply {$ixform} \
    -prefix {$T1_al}          \

# mask uncorr T1
set T1_mask = {$useDir}T1_{$analysis}_mask.nii.gz
3dAutomask -overwrite  \
    -prefix {$T1_mask} \
    {$T1_al}

# mask uncorr epi
set epi_mask = {$useDir}epi_{$analysis}_mask.nii.gz
3dAutomask -overwrite -prefix {$epi_mask} \
    {$use_epi}