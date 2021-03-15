#!/bin/csh
set whichAnalysis = 'HCP'
set topDir = '/*****PATH TO ANALYSIS DIRECTORY*****/'{$whichAnalysis}'/' # n.b. this must be set before running!
set scriptDir = {$topDir}'scripts/'

set subjNums = ('100610'\
    '102311'\
    '102816'\
    '104416'\
    '105923'\
    '108323'\
    '109123'\
    '111312'\
    '111514'\
    '114823'\
    '115017'\
    '115825'\
    '116726'\
    '118225'\
    '125525'\
    '126426'\
    '126931'\
    '128935'\
    '130114'\
    '130518')

@ idx = 1
while ( $idx <= $#subjNums )
    set subj = $subjNums[$idx]
    tcsh {$scriptDir}/do_analysis.sh $subj
    @ idx += 1
end