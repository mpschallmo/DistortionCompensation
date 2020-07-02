#!/bin/csh
set whichAnalysis = 'SBRef'
set topDir = '/*****PATH TO ANALYSIS DIRECTORY*****/'{$whichAnalysis}'/' # n.b. this must be set before running!
set scriptDir = {$topDir}'scripts/'

set subjNums = ('P6003691'\
    'P1010228'\
    'P6001501'\
    'P5104604'\
    'P6004604'\
    'P6004202'\
    'P1010299'\
    'P6004687'\
    'P1007451'\
    'P2104777'\
    'P6010671'\
    'P4100631'\
    'P2110465'\
    'P1010422'\
    'P6010465'\
    'P6004777'\
    'P1010407'\
    'P6010731'\
    'P6010363'\
    'P1006397'\
    'P4110363'\
    'P4104604'\
    'P3110692'\
    'P6010932'\
    'P3102476'\
    'P3111176'\
    'P1011139'\
    'P6004002'\
    'P1011033'\
    'P1010859'\
    'P1011399')
set scanTypes = ('Z'\
    'B'\
    'B'\
    'B'\
    'Z'\
    'Z'\
    'B'\
    'B'\
    'B'\
    'B'\
    'B'\
    'B'\
    'B'\
    'B'\
    'Z'\
    'Z'\
    'B'\
    'Z'\
    'Z'\
    'B'\
    'B'\
    'B'\
    'B'\
    'B'\
    'B'\
    'B'\
    'B'\
    'Z'\
    'B'\
    'B'\
    'B')

@ idx = 1
while ( $idx <= $#subjNums )
    set subj = $subjNums[$idx]
    set scan = $scanTypes[$idx]
    tcsh {$scriptDir}/do_one_scan.sh $subj $scan $whichAnalysis
    @ idx += 1
end