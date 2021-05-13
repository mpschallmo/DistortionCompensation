if ( $#argv == 3 ) then
    set subjDir = $argv[1]
    set analysis = $argv[2]
    set epi_num = $argv[3]
else
    echo "Please specify subjDir, analysis, and epi_num as an input variables"
endif

# this assumes you're sitting in a directory with all the files you need:
set rawDir = {$subjDir}'/raw/'
set outDir = {$subjDir}'/'{$analysis}'/head_motion/GE_'{$epi_num}'/'
set epi = 'epi_'{$analysis}'_'{$epi_num}
set ext = '.nii.gz'
set anat = '3T_anat_uni'
set anat_seg = {$rawDir}/'wmparc'
set dicePath = {$outDir}'/DICE.txt'

cd {$outDir}

rm $dicePath
touch $dicePath
echo "Dice coefficients:" >> $dicePath
echo "" >> $dicePath

set wr = {$outDir}/wmparc_{$analysis}_{$epi_num}
set xform = {$outDir}/epi_{$analysis}_{$epi_num}_al_3T_mat.aff12.1D
set ixform = {$outDir}/inv_epi_{$analysis}_{$epi_num}_al_3T_mat.aff12.1D

# get rid of any negative numbers (topup makes them)
3dcalc -overwrite -prefix {$epi}{$ext} -a {$epi}{$ext} -expr 'step(a)*a'

rm $ixform
cat_matvec -ONELINE {$xform} -I > $ixform

set anatFile = {$outDir}/{$anat}{$ext}
set anat_al = {$outDir}/{$anat}_al_EPI{$ext}

# mps 20191127 now warp 3T anat to match 7T EPI
3dAllineate                   \
    -overwrite                \
    -source {$anatFile}       \
    -base {$epi}{$ext}        \
    -1Dmatrix_apply {$ixform} \
    -prefix {$anat_al}        \

set anat_mask = {$outDir}/{$anat}_al_EPI_mask{$ext}
3dAutomask -overwrite    \
    -prefix {$anat_mask} \
    {$anat_al}
# mps 20191127 adding mask for 3T anat

# mps 20200108 mask epi before, so we can use for 3dSeg
3dAutomask -overwrite -prefix {$epi}_mask{$ext} \
    {$epi}{$ext}
    # mps 20191128 mask the dist comp epi

# mps 20191228 make wmparc_al manually instead of in minn.py
3dAllineate -overwrite        \
    -prefix wmparc_al{$ext}   \
    -1Dmatrix_apply {$ixform} \
    -master {$epi}{$ext}      \
    -final NN                 \
    {$anat_seg}{$ext}         \

# use wmparc_al from that to make GM, WM and CSF masks in EPI space
# labels from https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT
3dcalc -overwrite -prefix {$wr}_gm{$ext} -a wmparc_al.nii.gz    \
    -expr 'equals(a, 3) + equals(a, 8) + equals(a, 9) +         \
    equals(a, 10) + equals(a, 11) + equals(a, 12) +             \
    equals(a, 13) + equals(a, 17) + equals(a, 18) +             \
    equals(a, 19) + equals(a, 20) + equals(a, 26) +             \
    equals(a, 27) + equals(a, 30) + equals(a, 42) +             \
    step(a - 46)*step(57 - a) + equals(a, 58) +                 \
    equals(a, 59) + equals(a, 62) + equals(a, 96) +             \
    equals(a, 97) + step(a - 100)*step(108 - a) +               \
    step(a - 109)*step(117 - a) + step(a - 135)*step(140 - a) + \
    step(a - 192)*step(213 - a) + step(a - 213)*step(219 - a) + \
    equals(a, 220) + equals(a, 222) +                           \
    step(a - 399)*step(440 - a) + step(a - 499)*step(508 - a) + \
    step(a - 549)*step(558 - a) + step(a - 600)*step(629 - a) + \
    step(a- 639)*step(692-a) + step(a- 999)*step(3000-a)'

3dcalc -overwrite -prefix {$wr}_wm{$ext} -a wmparc_al.nii.gz   \
    -expr 'equals(a, 2) + equals(a, 7) +                       \
    equals(a, 16) + equals(a, 28) +                            \
    equals(a, 41) + equals(a, 46) + equals(a, 60) +            \
    equals(a, 77) + equals(a, 78) + equals(a, 79) +            \
    equals(a, 100) + equals(a, 108) + equals(a, 109) +         \
    equals(a, 117) + equals(a, 156) + equals(a, 157) +         \
    equals(a, 158) + equals(a, 192) +                          \
    step(a - 169)*step(176 - a) + equals(a, 219) +             \
    equals(a, 223) + step(a - 249)*step(256 - a) +             \
    equals(a, 508) + equals(a, 558) + step(a-2999)*step(5003-a)'

3dcalc -overwrite -prefix {$wr}_gmwm{$ext} -a {$wr}_gm{$ext} \
    -b {$wr}_wm{$ext} -expr 'step(a + b)' # mps 20191228 moving this above csf to use for masking

3dmerge -overwrite -prefix temp_gmwm_blur{$ext} -1blur_fwhm 5 {$wr}_gmwm{$ext}

3dcalc -overwrite -prefix temp_gmwm_blur_mask{$ext} -a temp_gmwm_blur{$ext} -expr 'step(a-.2)'

3dcalc -overwrite -prefix temp_csf.nii.gz -a {$wr}_gmwm{$ext} -b temp_gmwm_blur_mask{$ext} -expr 'b*step(1-a)'

3dcalc -overwrite -prefix {$wr}_csf{$ext} -a wmparc_al.nii.gz -b temp_csf.nii.gz \
    -expr 'step( equals(a, 1) + equals(a, 4) + equals(a, 5)                      \
    + equals(a, 6) + equals(a, 14) + equals(a, 15)                               \
    + equals(a, 24) + equals(a, 31) + equals(a, 40)                              \
    + equals(a, 43) + equals(a, 44) + equals(a, 45)                              \
    + equals(a, 63) + equals(a, 72) + equals(a, 122)                             \
    + equals(a, 213) + equals(a, 221) + equals(a, 256)                           \
    + step(b) )'
    # mps 20191128 changing the labels to be more thorough, include cerebellum and
    # subcortical structures
    # mps 20191228 changing CSF to include blurred edge from GMWM, as in segmentEPI.py
    # the blur gets the edge of the brain, and the parts from wmparc get bits in the ventricles

3dcalc -overwrite -prefix {$wr}_cset{$ext} -a {$wr}_gm{$ext} -b {$wr}_wm{$ext} \
    -c {$wr}_csf{$ext} -expr '1*c + 2*a + 3*b' # mps 20191228 adding cset

mv wmparc_al.nii.gz {$outDir}/wmparc_al.nii.gz # mps added 20191017 clean up

3dUnifize -overwrite -prefix {$epi}_uni{$ext} -input {$epi}{$ext} # mps 20191228 adding unifize EPI

3dSeg -classes 'CSF ; GM ; WM' -bias_classes 'GM ; WM' \
    -overwrite                                         \
    -mask {$epi}_mask{$ext}                            \
    -cset {$wr}_cset{$ext}                             \
    -anat {$epi}_uni{$ext} 
    # mps 20191228 run 3dSeg manually, using the wmparc files from above
    # mps 20200108 use mask from un-unifized EPI, looks better

3dAFNItoNIFTI -overwrite -prefix {$epi}_seg{$ext} Segsy/Classes+orig # mps 20191228

3dMedianFilter -overwrite -prefix {$epi}_seg{$ext} -irad 1 {$epi}_seg{$ext} # mps 20191228

rm -rf Segsy
rm temp_*.nii # mps added 20191017 clean up
rm temp_*{$ext} # mps added 20191228 clean up

# and make matching volumes for epi
3dcalc -overwrite -prefix {$epi}_csf{$ext} \
    -a {$epi}_seg{$ext} -expr 'equals(a, 1)'
3dcalc -overwrite -prefix {$epi}_gm{$ext} \
    -a {$epi}_seg{$ext} -expr 'equals(a, 2)'
3dcalc -overwrite -prefix {$epi}_wm{$ext} \
    -a {$epi}_seg{$ext} -expr 'equals(a, 3)'
3dcalc -overwrite -prefix {$epi}_gmwm{$ext}  \
    -a {$epi}_gm{$ext} -b {$epi}_wm{$ext} \
    -expr 'step(a + b)'

echo {$analysis}":" >> $dicePath
3ddot -dodice {$anat_mask} {$epi}_mask.nii.gz >> $dicePath
3ddot -dodice {$wr}_gmwm.nii.gz {$epi}_gmwm.nii.gz >> $dicePath
3ddot -dodice {$wr}_gm.nii.gz {$epi}_gm.nii.gz >> $dicePath
3ddot -dodice {$wr}_csf.nii.gz {$epi}_csf.nii.gz >> $dicePath
echo "" >> $dicePath
