# CMakeLists for studio project with specific directory tree
include(helpers.cmake)
cmake_minimum_required(VERSION 2.8.9)
project(studio_generated)
set(CMAKE_CXX_FLAGS "-std=c++14")

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-undefined,error")
else()
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "$ENV{STUDIO_PROJECT_BUILD_PATH}/lib")

# reading in cmake cache variables from build script
set(in_src ${IS_IN_SOURCE})
set(work_dir ${WORKING_DIR})
set(build_dir_list ${BUILD_DIRS})

# include directories
find_inc_dirs()
foreach(inc_dir ${inc_dirs_list})
  include_directories(${inc_dir})
endforeach()

find_lib_dirs()
foreach(lib_prod ${lib_prod_list})
  file(GLOB so_list ${map_${lib_prod}}/*.so ${map_${lib_prod}}/*.dylib)
  foreach(so ${so_list})
    string(REGEX REPLACE "^.*/lib([^/]*)\\.(so|dylib)$" "\\1" var ${so})
    set(${var} ${so})
  endforeach()
endforeach()

# root special case!
file(GLOB root_lib_list $ENV{ROOTSYS}/lib/*.so $ENV{ROOTSYS}/lib/*.dylib)
foreach(root_so ${root_lib_list})
  string(REGEX REPLACE "^/.*lib(.*)\\.(so|dylib)$" "\\1" root_var ${root_so})
  set(root_${root_var} ${root_so})
endforeach()

# define the set of libraries needed to link a service or module
# library lib lists per plugin
set(service_lib_list)
set(module_lib_list)
set(source_lib_list)
include(libsForPlugins.cmake)
libsForPlugins()

# adding sources

# if passed an arg
if(NOT build_dir_list STREQUAL "no args")
  message(STATUS "specified dirs are: ${build_dir_list}")
  # add each build dir
  foreach(build_dir ${build_dir_list})
    add_subdirectory($ENV{STUDIO_PROJECT_SRC_PATH]/${build_dir})
  endforeach(build_dir ${build_dir_list})
  # if no args passed
elseif(build_dir_list STREQUAL "no args")
  message(STATUS "no specified dirs...")
  # if outside srcs
  if(NOT in_src)
    message(STATUS "not in src dir... building everything...")
    add_subdirectory($ENV{STUDIO_PROJECT_SRC_PATH})
  elseif(in_src)
    message(STATUS "somewhere in src... building current dir...")
    add_subdirectory(${work_dir})
  endif()
endif(NOT build_dir_list STREQUAL "no args")
