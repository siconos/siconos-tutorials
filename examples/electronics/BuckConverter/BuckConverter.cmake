find_package(BLASDEV)
target_link_libraries(${_EXE} PUBLIC BLAS::BLAS)
