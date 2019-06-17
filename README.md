# unBOXER
Scripts for expanding BOXER-formatted compressed covariance matices to fully-specified covarianvce matrix form.

Based upon the Ye Olde Fortranne algorithm presented in Section 11.8 of the [NJOY2016 manual](https://github.com/njoy/NJOY2016-manual/raw/master/njoy16.pdf#page=372), roughly translated into modern Python. (See the [NJOY project](https://github.com/njoy/) for more details.)

Currently takes in reactions to be extracted from BOXER format via a (hard-coded) input file (fetchCov.dat), with the syntax:

    idType matNum1 MT1 matNum2 MT2

where `matNum1` and `MT1` are the ENDF/B-formatted material number and reaction number, respectively. Note that `iType = 3` and `iType = 4` correspond to reaction covariances; i.e., `matNum2` and `MT2` must be specified in addition to `matNum1` and `MT1`

### TODO ###

 - Add in grouping information by MAT (i.e., # reactions per material ID before each set of covariance matrices)
 - Add in translation from MAT & MT to symbolic representations (e.g., `2225 103 -> Ti-46(n,p)Sc-46`)
 - Develop validation tests against known covariance data
