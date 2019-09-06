## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   ENABLE_TESTING()
##   INCLUDE(CTest)

set(CTEST_PROJECT_NAME "Siconos examples")
set(CTEST_NIGHTLY_START_TIME "20:00:00 CET")

# https is needed on cdash server side
set(CTEST_DROP_METHOD "https")
set(CTEST_DROP_SITE "my.cdash.org/")
# set(CTEST_DROP_SITE "cdash-tripop.inrialpes.fr")
set(CTEST_DROP_LOCATION "/submit.php?project=Siconos")
set(CTEST_DROP_SITE_CDASH TRUE)
