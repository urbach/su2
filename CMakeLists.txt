# Copyright © 2022 Carsten Urbach <urbach@hiskp.uni-bonn.de>
#             2022 Simone Romiti <simone.romiti@uni-bonn.de>
# Licensed under the GNU Public License version 3

cmake_minimum_required(VERSION 3.16.3)
project(su2 C CXX)

set (CMAKE_CXX_STANDARD 17) # c++ standard

## here we generate version.hh
## to track the git version of the executables
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
else()
  set(GIT_BRANCH "")
  set(GIT_COMMIT_HASH "")
endif()

message(STATUS "Git current branch: ${GIT_BRANCH}")
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")

message(STATUS "Generating version.hh")
configure_file(
  ${CMAKE_SOURCE_DIR}/version.hh.in
  ${CMAKE_BINARY_DIR}/generated/version.hh
  )
  
  include_directories(${CMAKE_SOURCE_DIR}/)
  include_directories(${CMAKE_SOURCE_DIR}/include/)
  include_directories(${CMAKE_BINARY_DIR}/generated/)
  
  
  # openMP threads parallelization
  find_package(OpenMP)
  
  
  # setting C and CXX compilation flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -D_USE_OMP_")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -D_USE_OMP_")
  
  find_package(yaml-cpp)
  include_directories(${YAML_CPP_INCLUDE_DIR} ${YAML_CPP_INCLUDE_DIR}/yaml-cpp/)

  find_package(xtl)
  include_directories(${xtl_INCLUDE_DIRS} ${xtl_INCLUDE_DIRS}/xtl/)

  find_package(xtensor)
  include_directories(${xtensor_INCLUDE_DIRS} ${xtensor_INCLUDE_DIRS}/xtensor/)


  if(${HIGHFIVE_FORMAT})
    find_package(HighFive REQUIRED)
  endif()
  
# We could just have two `add_executable` blocks listing almost all source
# files. However, CMake would then compile the source files twice, once for
# each executable. The rationale is that different flags can be used for
# different executables. Here we want to re-use this. One can either track the
# generated objects files or create a simple static library. The latter is done
# here.
set(${CMAKE_INSTALL_LIBDIR} ${CMAKE_INSTALL_PREFIX}/lib/)
add_library(${PROJECT_NAME}
    exp_gauge.cc
    print_program_options.cc
    parse_commandline.cc
    parse_input_file.cc
)

# Create the main program.
find_package(Boost COMPONENTS program_options filesystem system)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

# executables only in the installation directory
set(${CMAKE_BINARY_DIR} ${CMAKE_INSTALL_PREFIX}/bin/)
add_executable(u1-main main-u1.cpp)
add_executable(su2-main main-su2.cpp)

add_executable(su2-kramers kramers.cc)
add_executable(test-groups test.cc)
add_executable(scaling scaling.cc)
add_executable(try tryreduce.cc)

set(test_progs "links" "smearing" "yaml")

if(${HIGHFIVE_FORMAT})
  set(test_progs ${test_progs} "h5-gen" "h5-dump")
endif()

set(test_progs_fullname "")
foreach(target ${test_progs})
  add_executable(test-${target}.exe test/${target}.cpp)
  list(APPEND test_progs_fullname test-${target}.exe) # now FOO="a;b"
endforeach(target)


target_link_libraries(test-yaml.exe ${YAML_CPP_LIBRARIES})

foreach(target ${test_progs_fullname})
  if(${HIGHFIVE_FORMAT})
    target_link_libraries(${target} HighFive)
  endif()
endforeach(target)


target_link_directories(${CMAKE_PROJECT_NAME} PUBLIC ${YAML_CPP_LIBRARY_DIR})

foreach(target u1-main su2-main su2-kramers test-groups scaling try ${test_progs_fullname})
  target_link_libraries(${target} su2 ${YAML_CPP_LIBRARIES})

  if(Boost_FOUND)
    target_link_libraries(${target} su2 ${Boost_LIBRARIES})
  else()
    target_link_libraries(${target} su2 "-lboost_program_options -lboost_filesystem -lboost_system")
  endif()
endforeach(target)

# installing 

install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(TARGETS u1-main su2-main)
install(TARGETS su2-kramers test-groups scaling try)

# tests

foreach(target ${test_progs})
  install(TARGETS test-${target}.exe LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endforeach(target)

