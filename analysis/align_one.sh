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
set anat = '3T_anat'
set unifizeT1 = 1

cd {$outDir}

    # unifize T1
if ( $unifizeT1 == 1 ) then
  if  !(-f {$rawDir}/{$anat}_uni{$ext}) then
    3dUnifize -overwrite -prefix {$rawDir}/{$anat}_uni{$ext} {$rawDir}/{$anat}{$ext}
  endif
  set anat = {$anat}'_uni'
endif

align_epi_anat.py                          \
     -epi2anat                             \
     -anat {$rawDir}/{$anat}{$ext}         \
     -rigid_body                           \
     -suffix _al_3T                        \
     -epi {$outDir}/{$epi}{$ext}           \
     -epi_base 0                           \
     -epi_strip 3dAutomask                 \
     -ginormous_move                       \
     -volreg off                           \
     -tshift off                           \
     -anat_has_skull no                    \
     -overwrite
mv epi_*_al* {$outDir}/.
mv 3T_anat_uni_al_3T_e2a_only_mat.aff12.1D \
    {$outDir}/3T_anat_uni_al_{$analysis}_mat.aff12.1D
3dcopy -overwrite {$rawDir}/{$anat}{$ext} {$outDir}/{$anat}{$ext}