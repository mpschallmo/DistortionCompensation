# DistortionCompensation
Code for "Assessing methods for geometric distortion compensation in 7T gradient echo fMRI data"

Code authors: Michael-Paul Schallmo and Cheryl A. Olman

1. System requirements
- Requires MATLAB, AFNI, FSL, Python, gradunwarp, and FreeSurfer (as implemented in the HCP Workbench)
- And requires the following MATLAB toolboxes: Image Processing Toolbox, Statistics and Machine Learning Toolbox
- Also requires the following python modules: os, nibabel, numpy, subprocess
- Shell scripts written for and run in tcsh
- Tested using MATLAB version R2017b, AFNI version 18.2.04, FSL version 5.0.9, Python version 2.7, gradunwarp version 1.0.3, FreeSurfer version 5.3 (HCP workbench version 3.22.0), all run on Linux Red Hat version 6.10

2. Installation guide
- Unzip or clone the directory
- The following scripts must be edited manually based on your local data paths,
look for "# n.b. this must be set before running!":
    - separate_GE/copy_raw_data.sh
    - separate_GE/do_all_scans.sh
    - separate_GE/do_all_six_script.sh
    - separate_GE/do_one_scan.sh
    - single_GE/do_all_scans.sh
    - single_GE/do_all_six_script.sh
    - SBRef/copy_raw_data.sh
    - SBRef/do_all_scans.sh
    - SBRef/do_all_six_script.sh
    - SBRef/do_one_scan.sh
    - each of the 4 scripts in SBRef/copy_subj_scripts/
- Time to install: 5 min

3. Instructions for use
- To process the data from one of the three analyses, run the following batch script from the command line:
tcsh do_all_scans.sh
- This will call the script do_one_scan.sh, which in turn calls the scripts that carry out each analysis step
- There are 3 different versions of do_all_scans.sh, one in each of the analysis sub-folders, which will batch process the data for that particular analysis
Time to run = about 12 hours on our server

- To analyze the data, run the following in MATLAB (from inside the analysis sub-directory, or with it added to your MATLAB paths):
results = dist_comp_analysis( options )
- For a description of the options and results structures, please see the help string (type “help dist_comp_analysis” or check the top of the function)
Time to run = about 5 minutes if the mutual information data matrices are not saved to memory, < 30 sec if they have been saved previously

4. Notes
This code uses the following publicly available MATLAB functions (included in the repository):
plotSpread - mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot
mutInfo - mathworks.com/matlabcentral/fileexchange/35625-information-theory-toolbox
Credit for these functions goes to their original authors (not us).

This code is made available under the following license:
Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
https://creativecommons.org/licenses/by-nc/4.0/
