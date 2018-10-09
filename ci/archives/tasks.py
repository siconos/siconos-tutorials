# --- Config to run siconos examples ---

# Note FP/MB : this should be the only task(s) that run examples !!!
#        We may add later some 'examples-with-bullet' or 'examples-with-mechanisms' ... tasks later.
#        All tasks 'examples' should:
#         - conf, make and make install of siconos components (no tests!)
#         - conf, make and make test of siconos examples

# Case1 : siconos 'basics' components, numerics, kernel, control and related examples
siconos_light_examples = minimal_with_python.copy()(
    ci_config='examples_light',
    targets={'.': ['docker-build', 'docker-cmake', 'docker-make',
                   'docker-make-install', 'docker-make-clean'],
             'examples': ['docker-build', 'docker-ctest', 'docker-make-clean']},
    add_srcs=['examples'])

# Case2 : siconos with mechanics components and bullet + related examples
siconos_all_examples = minimal_with_python.copy()(
    ci_config='examples_all',
    add_pkgs=['bullet', 'h5py'],
    targets={'.': ['docker-build', 'docker-cmake', 'docker-make',
                   'docker-make-install', 'docker-make-clean'],
             'examples': ['docker-build', 'docker-ctest', 'docker-make-clean']},
    add_srcs=['examples'])

siconos_test_deb = SiconosCiTask(
    ci_config='examples',
    distrib='ubuntu:16.04',
    pkgs=['siconos'],
    srcs=['examples'])

siconos_test_rpm = SiconosCiTask(
    ci_config='examples',
    distrib='fedora:latest',
    pkgs=['siconos'],
    srcs=['examples'])
