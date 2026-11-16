## Resubmission

This is a patch release fixing the test failures shown on the CRAN check page,
caused by numerical changes in lme4 >= 2.0 (improved GLMM covariance computation,
see https://github.com/lme4/lme4/issues/1003). Only test tolerances were changed;
package code is unchanged.

## Test environments
* local macOS, R 4.4.1, with lme4 2.1.0
* win-builder (devel)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Reverse dependencies
We checked the one strong reverse dependency (aniDom): no new problems.
