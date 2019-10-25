**This code was made for aa application challenge (7 days but not full-time) and was not completed after the deadline. Thus there are still some missing features and bugs.**

# Continuous Reaction-Diffusion

This code computes the dynamics of continuous reaction-diffusion problem.

WARNING: The random initialization is not ready yet!


## Requirements

This code needs the following libraries to be compiled:

* Boost (only tested with versions 1.65 and 1.70)
* CUDA (only tested with version 10.1)


## Installation

Run the following commands:

    git clone https://github.com/adrien-berchet/continuous_rd.git
    cd continuous_rd/build
    cmake .. # use -DCMAKE_INSTALL_PREFIX:PATH=<path> option to use custom installation path
    make all doc install


## Usage

Usage details are given with the following command:

	./continuous_rd.exe --help
