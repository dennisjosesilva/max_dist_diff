project(EDT)

cmake_minimum_required(VERSION 3.10)

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


# add_executable(dtmt_diff src/dtmt_diff.cpp)
# target_link_libraries(dtmt_diff
#     gft 
#     Threads::Threads
#     ZLIB::ZLIB
#     OpenMP::OpenMP_CXX)

# Exemplos com a morphotree
add_executable(max_dist_attr src/max_dist_attr.cpp)
target_link_libraries(max_dist_attr
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