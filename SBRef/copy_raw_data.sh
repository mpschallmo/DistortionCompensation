#!/bin/csh

# the user may specify a single subject to run with
if ( $#argv > 0 ) then
    set subj = $argv[1] # subj ID, e.g., P9999999
    set scanType = $argv[2] # B or Z?
    set whichAnalysis = $argv[3] # single_GE or separate_GE?
else
    echo "Please specify subject, scan letter, and which analysis (single_GE or separate_GE) as input variables"
endif

if ($scanType == A) then
  set scanNum = 2
else if ($scanType == B) then
  set scanNum = 2
else if ($scanType == Z) then
  set scanNum = 3
else
  echo "$scanType" is not a valid scanType
  exit
endif

# set up folders
set topDir = '/*****PATH TO ANALYSIS DIRECTORY*****/'{$whichAnalysis}'/' # n.b. this must be set before running!
set scanDir = {$topDir}{$subj}'/'{$scanType}'/'
set rawDir = {$scanDir}'raw/'
mkdir -p $rawDir

set dataDir = '/*****PATH TO NIFTI DATA DIRECTORY*****/' # n.b. this must be set before running!
set dicomDir = '/*****PATH TO DICOM DATA DIRECTORY*****/' # n.b. this must be set before running!
set nSlices = 85

# then convert the data
foreach scan (oppPE_AP PRF1_AP CSS_task3_AP COP_task1_AP FieldMap AP_SE REV_AP_SE t1_mpr_tra_iso_ND)
  if ( ($scan == 'AP_SE') || ($scan == 'REV_AP_SE') ) then
    set nTRs = 3
    to3d -prefix {$scan}.nii.gz -session {$rawDir} -time:zt {$nSlices} {$nTRs} 0 FROM_IMAGE {$dicomDir}????????-ST00?-PHCP{$subj}_V{$scanNum}_7{$scanType}/*BOLD_{$scan}/*.dcm

  else if ($scan == 'oppPE_AP') then
    3dcopy {$dataDir}'BOLDoppPEAP.nii.gz' {$rawDir}{$scan}.nii.gz

  else if ($scan == 'PRF1_AP') then
    3dcopy {$dataDir}'BOLDPRF1AP.nii.gz' {$rawDir}{$scan}.nii.gz
    set nTRs = 1
    to3d -prefix {$scan}_SBRef.nii.gz -session {$rawDir} -time:zt {$nSlices} {$nTRs} 0 FROM_IMAGE {$dicomDir}????????-ST00?-PHCP{$subj}_V{$scanNum}_7{$scanType}/*BOLD_{$scan}_SBRef/*.dcm

  else if ($scan == 'CSS_task3_AP') then
    3dcopy {$dataDir}'BOLDCSStask3AP.nii.gz' {$rawDir}{$scan}.nii.gz
    set nTRs = 1
    to3d -prefix {$scan}_SBRef.nii.gz -session {$rawDir} -time:zt {$nSlices} {$nTRs} 0 FROM_IMAGE {$dicomDir}????????-ST00?-PHCP{$subj}_V{$scanNum}_7{$scanType}/*BOLD_{$scan}_SBRef/*.dcm

  else if ($scan == 'COP_task1_AP') then
    3dcopy {$dataDir}'BOLDCOPtask1AP.nii.gz' {$rawDir}{$scan}.nii.gz
    set nTRs = 1
    to3d -prefix {$scan}_SBRef.nii.gz -session {$rawDir} -time:zt {$nSlices} {$nTRs} 0 FROM_IMAGE {$dicomDir}????????-ST00?-PHCP{$subj}_V{$scanNum}_7{$scanType}/*BOLD_{$scan}_SBRef/*.dcm
    
  else if ($scan == 'FieldMap') then
    set nTRs = 2
    to3d -prefix b0_mag.nii.gz -session {$rawDir} {$dicomDir}????????-ST00?-PHCP{$subj}_V{$scanNum}_7{$scanType}/MR*{$scan}/*-EC*.dcm
    to3d -prefix b0_ph.nii.gz -session {$rawDir} {$dicomDir}????????-ST00?-PHCP{$subj}_V{$scanNum}_7{$scanType}/PH*{$scan}/*-EC*.dcm

  else if ($scan == 't1_mpr_tra_iso_ND') then
    to3d -prefix {$scan}.nii.gz -session {$rawDir} {$dicomDir}????????-ST00?-PHCP{$subj}_V{$scanNum}_7{$scanType}/*{$scan}/*.dcm

  endif
end

# also grab 3T anatomy
set dir3T = '/*****PATH TO FREESURFER DATA DIRECTORY*****/PHCP'{$subj}'/T1w/PHCP'{$subj}'/mri/' # n.b. this must be set before running!
mri_convert {$dir3T}/brainmask.mgz {$rawDir}/3T_anat.nii.gz
mri_convert {$dir3T}/wmparc.mgz {$rawDir}/wmparc.nii.gz