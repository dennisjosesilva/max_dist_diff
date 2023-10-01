cmake_minimum_required(VERSION 3.10)

project(EDT)

set(CMAKE_CXX_STANDARD 14)

if(APPLE)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
    endif()

endif()

include_directories(
  "${CMAKE_CURRENT_SOURCE_DIR}/external/gft/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/external/morphotree/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/external/stb"
  "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/gft")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/external/morphotree")

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# Exemplos do professor Paulo
add_executable(dtmt src/dtmt.cpp)
target_link_libraries(dtmt
  gft 
  Threads::Threads
  ZLIB::ZLIB
  OpenMP::OpenMP_CXX)


add_executable(dtmt_diff src/dtmt_diff.cpp)
target_link_libraries(dtmt_diff
    gft 
    Threads::Threads
    ZLIB::ZLIB
    OpenMP::OpenMP_CXX)

# Exemplos com a morphotree
add_executable(max_dist_extiction_values 
  src/max_dist_extiction_values.cpp
  src/MaxDistComputer.cpp)

target_link_libraries(max_dist_extiction_values
  morphotree
  gft 
  Threads::Threads
  ZLIB::ZLIB
  OpenMP::OpenMP_CXX)

# display maximum distance attribute.
add_executable(max_dist_attr_disp 
  src/MaxDistComputer.cpp
  src/max_dist_attr_disp.cpp)
  
target_link_libraries(max_dist_attr_disp
  morphotree
  gft
  Threads::Threads
  ZLIB::ZLIB
  OpenMP::OpenMP_CXX)

add_executable(ultimate_maxDist_open 
src/MaxDistComputer.cpp
src/ultimate_maxDist_open.cpp)

target_link_libraries(ultimate_maxDist_open
morphotree
gft 
Threads::Threads
ZLIB::ZLIB
OpenMP::OpenMP_CXX)


# Direct filter
add_executable(direct_filter 
  src/BasicAttributeComputer.cpp
  src/direct_filter.cpp)
  
target_link_libraries(direct_filter
  morphotree  
  OpenMP::OpenMP_CXX)

# basic attr display
add_executable(basic_attr_disp 
  src/BasicAttributeComputer.cpp
  src/basic_attr_disp.cpp)
  
target_link_libraries(basic_attr_disp
  morphotree  
  OpenMP::OpenMP_CXX)  