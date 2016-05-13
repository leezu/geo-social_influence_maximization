# Geo-social influence maximization
Code for the bachelor thesis on geo-social influence maximization


## Coding
Code is written using features of C++14 standard.

It uses the following libraries:
- Boost graph library to store the graph and to facilitate graph algorithms
	- http://www.boost.org/doc/libs/1_60_0/libs/graph
	- Tutorial: https://github.com/richelbilderbeek/BoostGraphTutorial
- Boost spirit X3 to facilitate file parsing
	- http://ciere.com/cppnow15/x3_docs/
- docopt to generate the commandline interface
	- https://github.com/docopt/docopt.cpp


## Development
Clone the code with

    git clone --recursive git@github.com:leezu/geo-social_influence_maximization.git

To build the software make sure that cmake and the boost graph library are installed.
Debian/ubuntu packages are:
- cmake
- libboost-graph-dev

Boost Spirit X3 and docopt are included in the `lib` folder of the git repo, so they don't need to be installed.

Then go to an empty directory and run `cmake path/to/source_code` and `make` to build the code.
It creates an executable `gsinfmax` in the current folder.
