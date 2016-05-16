cmake_minimum_required(VERSION 3.2)

# Export compile commands for YouCompleteMe
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# The above is supposed to be enough to force C++14,
# but doesn't work for me
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

project(geo-social_influence_maximization VERSION 0.1 LANGUAGES CXX)

# Include boost graph library
# Load settings for the external project "Boost".
# The component "graph" is required.
find_package (Boost REQUIRED COMPONENTS graph)
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

# Include boost spirit X3 parser
set(Spirit_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/lib/spirit/include)
include_directories(INTERFACE ${Spirit_INCLUDE_DIR})
set(Spirit_HEADERS ${Spirit_INCLUDE_DIR}/boost/spirit/home/x3.hpp)

# Include docopt command line parsing tool
add_subdirectory(lib/docopt.cpp)

# "Add an executable to the project using the specified source files."
add_executable(gsinfmax main.cpp)

# "Link a target to given libraries."
target_link_libraries(gsinfmax docopt ${Boost_LIBRARIES} ${Spirit_HEADERS})
add_dependencies(gsinfmax docopt)

install(TARGETS gsinfmax DESTINATION bin)
