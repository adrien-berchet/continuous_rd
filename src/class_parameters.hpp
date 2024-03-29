/**
 * \file class_parameters.hpp
 *
 * \brief Provides a class to parse and store command line arguments.
 *
*/
#ifndef _CLASS_PARAMETERS
#define _CLASS_PARAMETERS

#include <iostream>
#include <fstream>

#include "program_options.hpp"


/**
 * \brief Class to process and store the parameters.
*/
class Parameters:public generic::ProgramOptions
{
public:
	/**
	 * \brief Define the parameters.
	*/
	void define_parameters()
	{
		code_name() = "RD_challenge";
		description() = "Compute Reaction-Diffusion process.";
		author() = "Adrien Berchet <adrien.berchet@gmail.com>";
		synopsis() = code_name() + " -p [FILE] [OPTIONS]...";

		add_param("-cu", "The coefficient cu", cu, 1, true);
		add_param("-cv", "The coefficient cv", cv, 1, true);
		add_param("-cw", "The coefficient cw", cw, 1, true);

		add_param("-c1", "The coefficient c1", c1, 1, true);
		add_param("-c2", "The coefficient c2", c2, 1, true);
		add_param("-c3", "The coefficient c3", c3, 1, true);
		add_param("-c4", "The coefficient c4", c4, 1, true);
		add_param("-c5", "The coefficient c5", c5, 1, true);
		add_param("-c6", "The coefficient c6", c6, 1, true);
		add_param("-c7", "The coefficient c7", c7, 1, true);
		add_param("-c8", "The coefficient c8", c8, 1, true);
		add_param("-c9", "The coefficient c9", c9, 1, true);

		add_param("-Du", "The coefficient Du", Du, 1, true);
		add_param("-Dv", "The coefficient Dv", Dv, 1, true);
		add_param("-Dw", "The coefficient Dw", Dw, 1, true);

		add_param("-Fmax", "The coefficient Fmax", Fmax, 1, true);
		add_param("-Gmax", "The coefficient Gmax", Gmax, 1, true);
		add_param("-Hmax", "The coefficient Hmax", Hmax, 1, true);

		add_param("-coeff_P", "The coefficient P", P, 1, true);

		add_param("-Nx", "The grid size along X axis", Nx, 1, true);
		add_param("-Ny", "The grid size along X axis", Ny, 1, true);
		add_param("-epsilon", "The grid step", epsilon, 1, false);
		add_param("-S", "The hexagon side length (the lattice is not generated if this argument is not provided)", S, 1, false);

		add_param("-gauss_std", "The standard deviation used for gaussian initialization", gauss_std, 1, false);

		add_param("-tmax", "The maximum time", tmax, 1, true);
		add_param("-dt", "The time step", dt, 1, true);

		add_param("-delta_obs", "The time step between two consecutive exports", delta_obs, 1, false);
		add_param("-result_folder", "The folder in which the results are exported", result_folder, 1, true);
		add_param("-export_neighbors", "Export neighbors for validation", export_neighbors, 1, false);
		add_param("-export_hex_lattice", "Export hexagonal lattice for validation", export_hex_lattice, 1, false);

		add_param("-log_level", "Level of logger verbosity (must be one of 'trace', 'debug', 'info', 'warning', 'error', 'fatal' or 'quiet')", log_level, 1, false);
	};

	/**
	 * \brief Constructor. Initialize default values and define parameters.
	*/
	Parameters()
	{
		epsilon=1;
		S=20. / sqrt(3.0);
		gauss_std=-1;
		export_neighbors=false;
		export_hex_lattice=false;
		log_level="debug";
		define_parameters();
	};

	/**
	 * \brief Copy constructor.
	*/
	Parameters(const Parameters &a)
	{
		cu=a.cu;
		cv=a.cv;
		cw=a.cw;
		c1=a.c1;
		c2=a.c2;
		c3=a.c3;
		c4=a.c4;
		c5=a.c5;
		c6=a.c6;
		c7=a.c7;
		c8=a.c8;
		c9=a.c9;
		Du=a.Du;
		Dv=a.Dv;
		Dw=a.Dw;
		Fmax=a.Fmax;
		Gmax=a.Gmax;
		Hmax=a.Hmax;
		P=a.P;
		tmax=a.tmax;
		Nx=a.Nx;
		Ny=a.Ny;
		dt=a.dt;
		epsilon=a.epsilon;
		S=a.S;
		gauss_std=a.gauss_std;
		delta_obs=a.delta_obs;
		result_folder=a.result_folder;
		export_neighbors=a.export_neighbors;
		export_hex_lattice=a.export_hex_lattice;
		log_level=a.log_level;
	};

	/**
	 * \brief Destructor.
	*/
	~Parameters() {};

	/**
	 * \brief Method to write the parameters into a file.
	*/
	void write(const std::string &name, const std::ios::openmode mode=std::ios::trunc)
	{
		std::ofstream file;
		file.open(name.c_str(), mode);
		file<<*this;
		file.close();
	};

public:
	/**@{ */
	// Model parameters
	double cu, cv, cw;
	double c1, c2, c3, c4, c5, c6, c7, c8, c9;
	double Du, Dv, Dw;
	double Fmax, Gmax, Hmax;
	double P;
	double tmax;
	size_t Nx, Ny;
	double dt;
	double epsilon;
	double S;
	double gauss_std;

	// IO parameters
	double delta_obs;
	std::string result_folder;
	bool export_neighbors;
	bool export_hex_lattice;
	std::string log_level;
	/**@}*/
};

#endif
