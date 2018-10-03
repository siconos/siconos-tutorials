Continuous integration configuration
====================================


Steps :

1. Download, build and install siconos 

   Described in install_siconos.sh
   
   Keep directory install-siconos (artifacts) for examples
   
   
2. Conf, build and run examples using ctest, with a submission to cdash.
   
   Described in build_and_run_examples.sh
   
   Use siconos install from previous step.
   
   Ctest process is described in ci/ctest_driver.cmake.
