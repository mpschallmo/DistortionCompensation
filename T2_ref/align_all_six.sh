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
set anat = '3T_anat' # Here, this is the T2 anat, mps 20210112
set unifizeT1 = 1

    # unifize T1
if ( $unifizeT1 == 1 ) then
  3dUnifize -overwrite -prefix {$rawDir}/{$anat}_uni{$ext} {$rawDir}/{$anat}{$ext}
  set anat = {$anat}'_uni'
endif

foreach analysis ('uncorr' 'fugue' 'GE_topup' 'SE_topup' 'GE_qwarp' 'SE_qwarp')
  set outDir = {$subjDir}/{$analysis}

  align_epi_anat.py                          \
       -epi2anat                             \
       -anat {$rawDir}/{$anat}{$ext}         \
       -rigid_body                           \
       -suffix _al_3T                        \
       -epi {$subjDir}/epi_{$analysis}{$ext} \
       -epi_base 0                           \
       -epi_strip 3dAutomask                 \
       -ginormous_move                       \
       -volreg off                           \
       -tshift off                           \
       -anat_has_skull no                    \
       -cost lpa                             \
       -overwrite
       # using cost lpa to align EPI to T2 anat mps 20201012
  mv epi_*_al* {$outDir}/.
  mv 3T_anat_uni_al_3T_e2a_only_mat.aff12.1D \
      {$outDir}/3T_anat_uni_al_{$analysis}_mat.aff12.1D
  3dcopy -overwrite {$rawDir}/{$anat}{$ext} {$outDir}/{$anat}{$ext}
end


set analysis = '12param'
set outDir = {$subjDir}/{$analysis}
mkdir $outDir
# same as above, but without -rigid_body, so 12 param instead
align_epi_anat.py -epi2anat           \
     -anat {$rawDir}/{$anat}{$ext}    \
     -suffix _al_3T                   \
     -epi {$subjDir}/epi_uncorr{$ext} \
     -epi_base 0                      \
     -epi_strip 3dAutomask            \
     -ginormous_move                  \
     -volreg off                      \
     -tshift off                      \
     -anat_has_skull no               \
     -cost lpa                        \
     -overwrite
     # using cost lpa to align EPI to T2 anat mps 20201012
mv epi_uncorr_al_3T+orig.BRIK {$outDir}/epi_{$analysis}_al_3T+orig.BRIK
mv epi_uncorr_al_3T+orig.HEAD {$outDir}/epi_{$analysis}_al_3T+orig.HEAD
mv epi_uncorr_al_3T_mat.aff12.1D {$outDir}/epi_{$analysis}_al_3T_mat.aff12.1D
mv 3T_anat_uni_al_3T_e2a_only_mat.aff12.1D \
    {$outDir}/3T_anat_uni_al_{$analysis}_mat.aff12.1D
3dcopy -overwrite {$rawDir}/{$anat}{$ext} {$outDir}/{$anat}{$ext}
3dcopy -overwrite {$subjDir}/epi_uncorr{$ext} {$outDir}/epi_12param{$ext}
# copy the uncorrected EPI because that's the "base" for the 12 param solution