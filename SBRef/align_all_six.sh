if ( $#argv > 0 ) then
    set subjDir = $argv[1]
else
    echo "Please specify subjDir as an input variable"
endif

# this assumes you're sitting in a directory with all the files you need:
set rawDir = {$subjDir}'/raw/'
set epi = 'PRF1_AP'
set GErev = 'oppPE_AP'
set SEfwd = 'AP_SE'
set SErev = 'REV_AP_SE'
set ext = '.nii.gz'
set anat = '3T_anat'
set unifizeT1 = 1

    # unifize T1
if ( $unifizeT1 == 1 ) then
  3dUnifize -overwrite -prefix {$rawDir}/{$anat}_uni{$ext} {$rawDir}/{$anat}{$ext}
  set anat = {$anat}'_uni'
endif

foreach analysis ('uncorr' 'fugue' 'GE_topup' 'SE_topup' 'GE_qwarp' 'SE_qwarp')
  set outDir = {$subjDir}/{$analysis}

  align_epi_anat.py                                \
       -epi2anat                                   \
       -anat {$rawDir}/{$anat}{$ext}               \
       -rigid_body                                 \
       -suffix _al_3T                              \
       -epi {$outDir}/epi_{$analysis}_SBRef{$ext}  \
       -epi_base 0                                 \
       -epi_strip 3dAutomask                       \
       -ginormous_move                             \
       -volreg off                                 \
       -tshift off                                 \
       -anat_has_skull no                          \
       -overwrite

  3dAllineate -overwrite                                     \
    -prefix {$outDir}/epi_{$analysis}_al_3T{$ext}            \
    -1Dmatrix_apply epi_{$analysis}_SBRef_al_3T_mat.aff12.1D \
    -source {$subjDir}/epi_{$analysis}{$ext}                 \
    -base {$outDir}/epi_{$analysis}_SBRef{$ext}              \

  mv epi_*_al* {$outDir}/.
  mv 3T_anat_uni_al_3T_e2a_only_mat.aff12.1D \
      {$outDir}/3T_anat_uni_al_{$analysis}_mat.aff12.1D
  3dcopy -overwrite {$rawDir}/{$anat}{$ext} {$outDir}/{$anat}{$ext}
end


set analysis = '12param'
set outDir = {$subjDir}/{$analysis}
mkdir $outDir
3dcopy -overwrite {$subjDir}/uncorr/epi_uncorr_SBRef{$ext} {$outDir}/epi_{$analysis}_SBRef{$ext}
# put the SBRef into the 12 param folder

# same as above, but without -rigid_body, so 12 param instead
align_epi_anat.py -epi2anat                     \
     -anat {$rawDir}/{$anat}{$ext}              \
     -suffix _al_3T                             \
     -epi {$outDir}/epi_{$analysis}_SBRef{$ext} \
     -epi_base 0                                \
     -epi_strip 3dAutomask                      \
     -ginormous_move                            \
     -volreg off                                \
     -tshift off                                \
     -anat_has_skull no                         \
     -overwrite

3dAllineate -overwrite                                     \
  -prefix {$outDir}/epi_{$analysis}_al_3T{$ext}            \
  -1Dmatrix_apply epi_{$analysis}_SBRef_al_3T_mat.aff12.1D \
  -source {$subjDir}/epi_uncorr{$ext}                      \
  -base {$outDir}/epi_{$analysis}_SBRef{$ext}              \

mv epi_{$analysis}_SBRef_al_3T+orig.BRIK {$outDir}/.
mv epi_{$analysis}_SBRef_al_3T+orig.HEAD {$outDir}/.
mv epi_{$analysis}_SBRef_al_3T_mat.aff12.1D {$outDir}/.
mv 3T_anat_uni_al_3T_e2a_only_mat.aff12.1D \
    {$outDir}/3T_anat_uni_al_{$analysis}_mat.aff12.1D
3dcopy -overwrite {$rawDir}/{$anat}{$ext} {$outDir}/{$anat}{$ext}
3dcopy -overwrite {$subjDir}/epi_uncorr{$ext} {$outDir}/epi_12param{$ext}
# copy the uncorrected EPI because that's the "base" for the 12 param solution