# define input variables
if ( $#argv > 0 ) then
    set subj = $argv[1] # subj ID, e.g., P9999999
    set scanType = $argv[2] # B or Z?
    set whichAnalysis = $argv[3] # single_GE, separate_GE?, or SBRef?
else
    echo "Please specify subject, scan letter, and which analysis (single_GE, separate_GE, or SBRef) as input variables"
endif


# set folders
set topDir = '/*****PATH TO ANALYSIS DIRECTORY*****/'{$whichAnalysis}'/' # n.b. this must be set before running!
set scriptDir = {$topDir}'scripts/'
set subjScriptDir = {$scriptDir}'copy_subj_scripts/'
set scanDir = {$topDir}{$subj}'/'{$scanType}'/'
set outFile = {$scanDir}'/out.'{$subj}'.'{$scanType}'.FMSD.txt'

mkdir -p $scanDir
rm $outFile
touch $outFile

# convert data
if !( -d $subjScriptDir ) then
    mkdir -p $subjScriptDir
endif
set subjScript = {$subjScriptDir}{$subj}'_'{$scanType}'_'{$whichAnalysis}'.sh'
if !( -f $subjScript ) then
    cp {$scriptDir}copy_raw_data.sh {$subjScript}
endif
tcsh {$subjScript} $subj $scanType $whichAnalysis |& tee $outFile

# do all six (distortion compensation methods)
tcsh {$scriptDir}do_all_six_script.sh $scanDir |& tee -a $outFile

# align all six, also does 12 param alingment of uncorrected
tcsh {$scriptDir}align_all_six.sh $scanDir |& tee -a $outFile

# segment all six
tcsh {$scriptDir}segment_all_six.sh $scanDir |& tee -a $outFile