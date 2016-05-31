cmake_minimum_required(VERSION 3.2)

# Export compile commands for YouCompleteMe
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable Release profile (and therefore optimization flags)
set(CMAKE_BUILD_TYPE Release)

project(geo-social_influence_maximization VERSION 0.1 LANGUAGES CXX)

# Include boost graph library
# Load settings for the external project "Boost".
# The component "graph" is required.
find_package (Boost REQUIRED COMPONENTS graph iostreams)
if (Boost_FOUND)
  # From the offical documentation:
  # Add include directories to the build. [...] If the SYSTEM option is given,
  # the compiler will be told the directories are meant as system include
  # directories on some platforms (signalling this setting might achieve effects
  # such as the compiler skipping warnings [...])."
  include_directories (SYSTEM ${Boost_INCLUDE_DIR})

  # From the offical documentation:
  # "Specify directories in which the linker will look for libraries. [...] Note
  # that this command is rarely necessary. Library locations returned by
  # find_package() and find_library() are absolute paths. Pass these absolute
  # library file paths directly to the target_link_libraries() command. CMake
  # will ensure the linker finds them."
  #link_directories (${Boost_LIBRARY_DIRS})
endif ()

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# Include boost spirit X3 parser
set(Spirit_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/lib/spirit/include)
include_directories(INTERFACE ${Spirit_INCLUDE_DIR})

# Include docopt command line parsing tool
include_directories(lib/docopt.cpp)
add_library(docopt STATIC lib/docopt.cpp/docopt.cpp)

# Add an executable to the project using the specified source files
set(gsinfmax_SRC main.cpp Graph_reader.cpp misc.cpp algorithms-lazy_greedy.cpp)
add_executable(gsinfmax ${gsinfmax_SRC})

# Link against boost libraries.
target_link_libraries(gsinfmax docopt ${Boost_LIBRARIES})

install(TARGETS gsinfmax DESTINATION bin)
