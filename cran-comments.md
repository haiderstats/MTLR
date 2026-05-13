## Test environments
* Local: macOS arm64, R 4.6.0
* GitHub Actions:
  * windows-latest: release
  * macOS-latest: release
  * ubuntu-latest: release
* R-hub:
  * ubuntu-latest (R-release)
  * windows-latest (R-release)
  * Fedora Linux 42 / atlas (R-devel)
  * Ubuntu 22.04 / ubuntu-clang (R-devel)
* win-builder: oldrelease, release, and devel

## R CMD check results
There were no ERRORs or WARNINGs. There is 1 NOTE:

* New maintainer: email updated from hshaider@ualberta.ca to humzahaider0@gmail.com.

## Resubmission
This is a resubmission of a failed 0.2.2 submission from 2019 (confirmation 
email was missed). Changes in this version:

* Removed spurious LazyData field from DESCRIPTION.
* Updated Authors@R field to replace deprecated Author/Maintainer fields.
* Replaced deprecated ggplot2 functions: aes_string() -> aes(), size -> linewidth.
* Replaced deprecated survival::survConcordance() with survival::concordance().

## Downstream dependencies
We checked 1 reverse dependency (CPSM 1.4.0), comparing R CMD check results
across CRAN and dev versions of this package.

* No new problems were detected.
