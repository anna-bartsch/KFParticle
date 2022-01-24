
# Instrucions:
#
# 1. type "make all test" to be sure that the library is rebuild before test
# 2. to see the text output, run the executables from test/. folder manually
#

# compile with the same flags as KFParticle to be sure that includes are the same

get_target_property(KFFlags "KFParticle" COMPILE_FLAGS)
set(CMAKE_CXX_FLAGS ${KFFlags})

message(STATUS "unit test: read compilation flags from the KFParticle library: ${KFFlags}")
message(STATUS "unit test: FIXTARGET flag is set to ${FIXTARGET}")

# make sure that the KFParticle library will be linked from this build, 
# not from LD_LIBRARY_PATH

link_libraries( KFParticle "-Wl,--disable-new-dtags")

# the tests
add_executable(kfparticle-unit-test test.cxx)
add_test(FirstTest kfparticle-unit-test)
