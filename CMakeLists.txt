cmake_minimum_required(VERSION 3.16)
project(zinc VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)
FetchContent_Declare(eigen3 URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz)
FetchContent_Declare(pugixml URL https://github.com/zeux/pugixml/archive/refs/tags/v1.14.tar.gz)
FetchContent_Declare(CLI11 URL https://github.com/CLIUtils/CLI11/archive/refs/tags/v2.4.2.tar.gz)
FetchContent_Declare(fmt URL https://github.com/fmtlib/fmt/archive/refs/tags/10.2.1.tar.gz)
FetchContent_Declare(googletest URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz)

# Prevent dependency test suites from being added to this project's CTest run.
set(BUILD_TESTING OFF CACHE BOOL "Disable tests in fetched dependencies" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(gtest_build_tests OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(eigen3 pugixml CLI11 fmt googletest)

add_library(zinc_core
  src/core/constants.cpp
  src/core/lattice.cpp
  src/core/structure.cpp
  src/io/driver.cpp
  src/io/qe_xml.cpp
  src/io/qe_log.cpp
  src/io/qe_input.cpp
  src/io/poscar.cpp
  src/io/cif.cpp
  src/io/xyz.cpp
)
target_include_directories(zinc_core PUBLIC include)
target_link_libraries(zinc_core PUBLIC Eigen3::Eigen pugixml fmt::fmt)

add_executable(zinc_cli src/cli/main.cpp src/cli/commands/scan.cpp src/cli/commands/info.cpp src/cli/commands/convert.cpp src/cli/commands/update.cpp)
target_link_libraries(zinc_cli PRIVATE zinc_core CLI11::CLI11)
set_target_properties(zinc_cli PROPERTIES OUTPUT_NAME zinc)

enable_testing()
add_executable(zinc_tests tests/test_core.cpp tests/test_io.cpp tests/test_cli.cpp)
target_link_libraries(zinc_tests PRIVATE zinc_core GTest::gtest_main GTest::gmock)
add_test(NAME zinc_tests COMMAND zinc_tests)

install(TARGETS zinc_cli DESTINATION bin)
