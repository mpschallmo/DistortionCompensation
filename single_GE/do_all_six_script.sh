# the user may specify a single subject to run with
if ( $#argv > 0 ) then
    set subjDir = $argv[1]
else
    echo "Please specify subjDir as an input variable"
endif

# if you want to do a subset, 0 turns stuff off ...
set gradUNWARP = 1
set VOLREG = 1
set doFUGUE = 1
set GE_TOPUP = 1
set SE_TOPUP = 1
set b0_AFNI = 1
set GE_QWARP = 1
set SE_QWARP = 1
set doDice = 0
set uncorr = 1
set yes = 1

# this assumes you're sitting in a directory with all the files you need:
set rawDir = {$subjDir}'/raw'
set epi_forGE = 'PRF1_AP'
set epi_forb0 = $epi_forGE # only 1 EPI scan
set epi_forSE = $epi_forGE # only 1 EPI scan
set epi = $epi_forGE
set use_suffix = '' # this will be used to append processing notes to the end of file names
set GErev = 'oppPE_AP'
set SEfwd = 'AP_SE'
set SErev = 'REV_AP_SE'
set mag = 'b0_mag'
set ph = 'b0_ph'
set dwell = '0.00032' # this needs to be actual dwell / iPAT!
set TRO = '41.60'
set ext = '.nii.gz'
set coefFile = '/*****PATH TO GIT REPO*****/ref_files/7TAScoeff_fromScanner_12Jan2018.grad' # n.b. this must be set before running!
set topup_acq_params = '/*****PATH TO GIT REPO*****/ref_files/topup_acq_params.txt' # n.b. this must be set before running!
set scriptDir = '/*****PATH TO ANALYSIS DIRECTORY*****/'{$whichAnalysis}'/scripts' # n.b. this must be set before running!

# first pluck out 3 TRs of task scans, use closest to adjacent FM scans
3dTcat -overwrite -prefix {$rawDir}/{$epi}{$use_suffix}{$ext} {$rawDir}/{$epi}{$use_suffix}{$ext}'[0..2]'

# take the median
3dTstat -overwrite -median -prefix {$rawDir}/{$epi}{$use_suffix}{$ext} {$rawDir}/{$epi}{$use_suffix}{$ext}

# then collapse EPI FM scans to 1 TR median
3dTstat -overwrite -median -prefix {$rawDir}/{$GErev}{$use_suffix}{$ext} {$rawDir}/{$GErev}{$use_suffix}{$ext}
3dTstat -overwrite -median -prefix {$rawDir}/{$SEfwd}{$use_suffix}{$ext} {$rawDir}/{$SEfwd}{$use_suffix}{$ext}
3dTstat -overwrite -median -prefix {$rawDir}/{$SErev}{$use_suffix}{$ext} {$rawDir}/{$SErev}{$use_suffix}{$ext}
3dTstat -overwrite -median -prefix {$rawDir}/{$mag}{$use_suffix}{$ext} {$rawDir}/{$mag}{$use_suffix}{$ext}

# for the sake of thoroughness, apply grad_unwarp to all of them!
if ( $gradUNWARP == $yes) then

  set inFile = {$rawDir}/{$epi}{$use_suffix}{$ext}
  set outFile = {$rawDir}/rm_{$epi}{$use_suffix}_gu{$ext}
  gradient_unwarp.py $inFile $outFile siemens -g $coefFile -n
  3drefit -view orig -space ORIG $outFile
  python {$scriptDir}/makeFullWarpRel.py fullWarp_abs.nii.gz
  mv fullWarp_abs.nii.gz {$rawDir}/fullWarp_abs.nii.gz
  # this gets written in whatever the current working directory is so move it to raw
  mv fullWarp_rel.nii.gz {$rawDir}/fullWarp_rel.nii.gz

  set old_suffix = $use_suffix
  set use_suffix = $use_suffix'_gu'

  foreach scan ( $epi $GErev $SEfwd $SErev $mag $ph )
    set inFile = {$rawDir}/{$scan}{$old_suffix}{$ext}
    set use_warp = {$rawDir}/fullWarp_rel.nii.gz
    3dNwarpApply                                     \
      -source $inFile                                \
      -nwarp $use_warp \
      -prefix {$rawDir}/{$scan}{$use_suffix}{$ext}
  end
endif

# next align all 3 task scans, so we can align the FM scans...
if ( $VOLREG == $yes ) then
  set old_suffix = $use_suffix
  set use_suffix = $use_suffix'_al'
  foreach scan ( $epi $GErev $SEfwd $SErev )
    align_epi_anat.py                           \
      -overwrite                                \
      -rigid_body                               \
      -epi2anat                                 \
      -anat {$rawDir}/{$mag}{$old_suffix}{$ext} \
      -suffix _al                               \
      -epi {$rawDir}/{$scan}{$old_suffix}{$ext} \
      -epi_base 0                               \
      -dset1_strip 3dSkullStrip                 \
      -dset2_strip 3dAutomask                   \
      -big_move                                 \
      -volreg off                               \
      -tshift off                               \
      -cost lpa                                 \
      -anat_has_skull yes
  end

  mv *_al* {$rawDir}/.
endif

# now apply logic to resample in 1 step, if doing both gradunwarp and volreg
foreach scan ( $epi $GErev $SEfwd $SErev $mag $ph )
  if ( $gradUNWARP == $yes ) then
    set warpList = {$rawDir}/fullWarp_rel.nii.gz 
  else
    set warpList = ''
  endif

  if ( $VOLREG == $yes ) then
    if ( ( $scan == $mag ) || ( $scan == $ph ) ) then
      set warpList = {$warpList}
      # don't need to align to itself
    else
      set warpList = {$warpList}' '{$rawDir}/{$scan}{$use_suffix}_mat.aff12.1D
    endif
  endif

  set inFile = {$rawDir}/{$scan}{$ext}
  # explicitly NOT using suffix here, don't want to resample more than 1 time
  
  set outFile = {$rawDir}/{$scan}{$use_suffix}{$ext}
  # do want to apply suffixes here

  3dNwarpApply       \
    -overwrite       \
    -source $inFile  \
    -nwarp "$warpList" \
    -prefix $outFile

  endif
end

# apply suffix to file names
set epi_forGE = {$epi_forGE}{$use_suffix}
set epi_forb0 = {$epi_forb0}{$use_suffix}
set epi_forSE = {$epi_forSE}{$use_suffix}
set epi = {$epi}{$use_suffix}
set GErev = {$GErev}{$use_suffix}
set SEfwd = {$SEfwd}{$use_suffix}
set SErev = {$SErev}{$use_suffix}
set mag = {$mag}{$use_suffix}
set ph = {$ph}{$use_suffix}

# clean up rm files (temp)
rm -f {$rawDir}/rm*

# now we're ready to start distortion compensation

###################################################
## uncorrected
###################################################
if ( $uncorr == $yes ) then
  set outDir = {$subjDir}'/uncorr'
  mkdir $outDir
  3dcopy -overwrite {$rawDir}/{$epi}{$ext} {$outDir}/epi_uncorr{$ext}
  3dcopy -overwrite {$rawDir}/{$epi}{$ext} {$subjDir}/epi_uncorr{$ext}
endif

# ###################################################
# ## FUGUE
# ###################################################
if ($doFUGUE == $yes) then
  set outDir = {$subjDir}'/fugue'
  mkdir $outDir
  bet {$rawDir}/{$mag} {$outDir}/mag_bet -m 

  3dresample -overwrite -prefix {$outDir}/mag_bet{$ext} -input {$outDir}/mag_bet{$ext} \
    -master {$rawDir}/{$ph}{$ext}

  fsl_prepare_fieldmap SIEMENS {$rawDir}/{$ph}{$ext} {$outDir}/mag_bet{$ext} \
    {$outDir}/ph_rad_per_sec{$ext} 1.02
  3dMedianFilter -overwrite -prefix {$outDir}/ph_rad_per_sec_medfilt{$ext} \
    -irad 1 {$outDir}/ph_rad_per_sec{$ext}
  fugue -i {$rawDir}/$epi_forb0 --loadfmap={$outDir}/ph_rad_per_sec_medfilt{$ext} \
    --dwell={$dwell} --unwarpdir=y -u {$outDir}/epi_fugue

  cp {$outDir}/epi_fugue{$ext} {$subjDir}/.
endif

# ###################################################
# ## TOPUP - GE
# ###################################################
if ($GE_TOPUP == $yes) then
  set outDir = {$subjDir}'/GE_topup'
  mkdir $outDir
  set fwdScan = {$rawDir}/{$epi_forGE}{$ext}
  set revScan = {$rawDir}/{$GErev}{$ext}
  # create median datasets from forward and reverse time series
  set fwd_med = {$outDir}/fwd_med{$ext}
  set rev_med = {$outDir}/rev_med{$ext}
  3dTstat -overwrite -median -prefix $fwd_med $fwdScan
  3dTstat -overwrite -median -prefix $rev_med $revScan

  # topup needs extra slice
  set applyTo = {$rawDir}/{$epi_forGE}{$ext}
  set padded_epi = {$outDir}/padded_epi{$ext}
  3dZeropad -overwrite -prefix $fwd_med -S 1 $fwd_med
  3dZeropad -overwrite -prefix $rev_med -S 1 $rev_med
  3dZeropad -overwrite -prefix $padded_epi -S 1 $applyTo

  set both_vol = {$outDir}/both_volumes{$ext}
  3dTcat -overwrite -prefix $both_vol $fwd_med $rev_med

  set topup_map = {$outDir}/topup_map
  topup --imain=$both_vol --datain=$topup_acq_params --config=b02b0.cnf \
    --out=$topup_map
  applytopup --imain=$padded_epi --inindex=1 --datain=$topup_acq_params \
    --method=jac --topup=$topup_map --out={$outDir}/padded_topped_up
  3dZeropad -overwrite -prefix {$outDir}/epi_GE_topup{$ext} -S -1 \
    {$outDir}/padded_topped_up{$ext}
  cp {$outDir}/epi_GE_topup{$ext} {$subjDir}/.
endif

# ###################################################
# ## TOPUP - SE
# ###################################################
if ($SE_TOPUP == $yes) then
  set outDir = {$subjDir}'/SE_topup'
  mkdir $outDir
  set fwdScan = {$rawDir}/{$SEfwd}{$ext}
  set revScan = {$rawDir}/{$SErev}{$ext}
  # create median datasets from forward and reverse time series
  set fwd_med = {$outDir}/fwd_med{$ext}
  set rev_med = {$outDir}/rev_med{$ext}
  3dTstat -overwrite -median -prefix $fwd_med $fwdScan
  3dTstat -overwrite -median -prefix $rev_med $revScan

  # topup needs extra slice
  set applyTo = {$rawDir}/{$epi_forSE}{$ext}
  set padded_epi = {$outDir}/padded_epi{$ext}
  3dZeropad -overwrite -prefix $fwd_med -S 1 $fwd_med
  3dZeropad -overwrite -prefix $rev_med -S 1 $rev_med
  3dZeropad -overwrite -prefix $padded_epi -S 1 $applyTo

  set both_vol = {$outDir}/both_volumes{$ext}
  3dTcat -overwrite -prefix $both_vol $fwd_med $rev_med

  set topup_map = {$outDir}/topup_map
  topup --imain=$both_vol --datain=$topup_acq_params --config=b02b0.cnf \
    --out=$topup_map
  applytopup --imain=$padded_epi --inindex=1 --datain=$topup_acq_params \
    --method=jac --topup=$topup_map --out={$outDir}/padded_topped_up
  3dZeropad -overwrite -prefix {$outDir}/epi_SE_topup{$ext} -S -1 \
    {$outDir}/padded_topped_up{$ext}
  cp {$outDir}/epi_SE_topup{$ext} {$subjDir}/.
endif


# ###################################################
# ## QWARP - GE
# ###################################################
if ($GE_QWARP == $yes) then
  set outDir = {$subjDir}'/GE_qwarp'
  mkdir $outDir
  set fwdScan = {$rawDir}/{$epi_forGE}{$ext}
  set revScan = {$rawDir}/{$GErev}{$ext}
  # create median datasets from forward and reverse time series
  set fwd_med = {$outDir}/fwd_med{$ext}
  set rev_med = {$outDir}/rev_med{$ext}

  # create median datasets from forward and reverse time series
  3dTstat -median -prefix $fwd_med $fwdScan
  3dTstat -median -prefix $rev_med $revScan

  # automask the median datasets 
  3dAutomask -apply_prefix {$outDir}/fwd_masked.nii.gz $fwd_med
  3dAutomask -apply_prefix {$outDir}/rev_masked.nii.gz $rev_med

  # compute the midpoint warp between the median datasets
  3dQwarp -plusminus -pmNAMES Rev For                           \
        -pblur 0.05 0.05 -blur -1 -1                          \
        -noweight -minpatch 9                                 \
        -source {$outDir}/rev_masked.nii.gz                   \
        -base   {$outDir}/fwd_masked.nii.gz                   \
        -prefix {$outDir}/blip_warp

  3dNwarpApply -quintic -nwarp {$outDir}/blip_warp_For_WARP+orig      \
             -source {$rawDir}/{$epi_forGE}{$ext}         \
             -prefix {$outDir}/epi_GE_qwarp{$ext}
  3drefit -atrcopy {$rawDir}/{$epi_forGE}{$ext} IJK_TO_DICOM_REAL      \
                 {$outDir}/epi_GE_qwarp{$ext}
  cp {$outDir}/epi_GE_qwarp{$ext} {$subjDir}/.
endif

# ###################################################
# ## QWARP - SE
# ###################################################
if ($SE_QWARP == $yes) then
  set outDir = {$subjDir}'/SE_qwarp'
  mkdir $outDir
  set fwdScan = {$rawDir}/{$SEfwd}{$ext}
  set revScan = {$rawDir}/{$SErev}{$ext}
  # create median datasets from forward and reverse time series
  set fwd_med = {$outDir}/fwd_med{$ext}
  set rev_med = {$outDir}/rev_med{$ext}

  # create median datasets from forward and reverse time series
  3dTstat -median -prefix $fwd_med $fwdScan
  3dTstat -median -prefix $rev_med $revScan

  # automask the median datasets 
  3dAutomask -apply_prefix {$outDir}/fwd_masked.nii.gz $fwd_med
  3dAutomask -apply_prefix {$outDir}/rev_masked.nii.gz $rev_med

  # compute the midpoint warp between the median datasets
  3dQwarp -plusminus -pmNAMES Rev For                           \
        -pblur 0.05 0.05 -blur -1 -1                          \
        -noweight -minpatch 9                                 \
        -source {$outDir}/rev_masked.nii.gz                   \
        -base   {$outDir}/fwd_masked.nii.gz                   \
        -prefix {$outDir}/blip_warp

  3dNwarpApply -quintic -nwarp {$outDir}/blip_warp_For_WARP+orig      \
             -source {$rawDir}/{$epi_forSE}{$ext}         \
             -prefix {$outDir}/epi_SE_qwarp{$ext}
  3drefit -atrcopy {$rawDir}/{$epi_forSE}{$ext} IJK_TO_DICOM_REAL      \
                 {$outDir}/epi_SE_qwarp{$ext}
  cp {$outDir}/epi_SE_qwarp{$ext} {$subjDir}/.
endif

# Make composite images for easy viewing
set outName = 'input_series'
3dTcat -overwrite -prefix {$subjDir}/{$outName}{$ext} \
  {$rawDir}/{$epi_forGE}{$ext}                        \
  {$rawDir}/{$GErev}{$ext}                            \
  {$rawDir}/{$epi_forb0}{$ext}                        \
  {$rawDir}/{$mag}{$ext}                              \
  {$rawDir}/{$epi_forSE}{$ext}                        \
  {$rawDir}/{$SEfwd}{$ext}                            \
  {$rawDir}/{$SErev}{$ext}                            \

set outName = 'output_series'
3dTcat -overwrite -prefix {$subjDir}/{$outName}{$ext} \
  {$subjDir}/epi_uncorr{$ext}                         \
  {$subjDir}/epi_fugue{$ext}                          \
  {$subjDir}/epi_b0_afni{$ext}                        \
  {$subjDir}/epi_GE_topup{$ext}                       \
  {$subjDir}/epi_SE_topup{$ext}                       \
  {$subjDir}/epi_GE_qwarp{$ext}                       \
  {$subjDir}/epi_SE_qwarp{$ext}                        