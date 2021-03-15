# define input variables
if ( $#argv > 0 ) then
    set subj = $argv[1] # subj ID, e.g., P9999999
    set scanType = $argv[2] # B or Z?
    set whichAnalysis = $argv[3] # single_GE or separate_GE?
else
    echo "Please specify subject, scan letter, and which analysis (single_GE or separate_GE) as input variables"
endif


# set folders
set topDir = '/*****PATH TO ANALYSIS DIRECTORY*****/'{$whichAnalysis}'/' # n.b. this must be set before running!
set scriptDir = {$topDir}'scripts/'
set scanDir = {$topDir}{$subj}'/'{$scanType}'/'
set outFile = {$scanDir}'/out.'{$subj}'.'{$scanType}'.FMSD.txt'

mkdir -p $scanDir
rm $outFile
touch $outFile

# convert data
tcsh {$scriptDir}copy_raw_data.sh $subj $scanType $whichAnalysis |& tee $outFile

# do all six (distortion compensation methods)
tcsh {$scriptDir}do_all_six_script.sh $scanDir |& tee -a $outFile

# align all six, also does 12 param alingment of uncorrected
tcsh {$scriptDir}align_all_six.sh $scanDir |& tee -a $outFile

# segment all six
tcsh {$scriptDir}segment_all_six.sh $scanDir |& tee -a $outFile