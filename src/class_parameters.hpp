#ifndef _CLASS_PARAMETERS
#define _CLASS_PARAMETERS

#include <iostream>
#include <fstream>

#include "program_options.hpp"

static const long double pi=3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117L;
static const long double piinv=1.L/pi;
static const long double twopi=2.L*pi;
static const long double twopi_inv=1.L/twopi;


// **********************************************************************
//		classe contenant tous les paramètres du transport
// **********************************************************************
class Parameters:public generic::ProgramOptions
{
public:
	void define_parameters()
	{
		code_name() = "RD_challenge";
		description() = "Compute Reaction-Diffusion process.";
		author() = "Adrien Berchet <adrien.berchet@gmail.com>";
		synopsis() = code_name() + " -i [FILE] -o [FILE] [OPTIONS]...";

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
		add_param("-tmax", "The maximum time", tmax, 1, true);
		add_param("-Nx", "The grid size along X axis", Nx, 1, true);
		add_param("-Ny", "The grid size along X axis", Ny, 1, true);
		add_param("-dt", "The time step", dt, 1, true);
		add_param("-epsilon", "The coefficient epsilon", epsilon, 1, true);
		add_param("-S", "The hexagon side length (the lattice is not generated if this argument is not provided)", S, 1, false);
		add_param("-delta_obs", "The time step between two consecutive exports", delta_obs, 1, true);
		add_param("-result_folder", "The folder in which the results are exported", result_folder, 1, true);
	};

	Parameters()
	{
		define_parameters();
	};

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
		delta_obs=a.delta_obs;
		result_folder=a.result_folder;
	};

	~Parameters() {};

	void write(const std::string &name, const std::ios::openmode mode=std::ios::trunc)
	{
		std::ofstream file;
		file.open(name.c_str(), mode);
		file<<*this;
		file.close();
	};

public:
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

	// IO parameters
	double delta_obs;
	std::string result_folder;
};

// template <typename T>
// std::ostream& operator<<(std::ostream &out, const params<T> &a)
// {
// 	out<<std::endl;
// 	out<<"********************************************"<<std::endl;
// 	out<<"*************** Parameters *****************"<<std::endl;
// 	out<<"********************************************"<<std::endl;
// 	out<<"cu="<<cu<<"\tcv="<<cv<<"\tcw="<<cw<<std::endl;
//     out<<"c1="<<c1<<"\tc2="<<c2<<"\tc3="<<c3<<"\tc4="<<c4<<"\tc5="<<c5<<"\tc6="<<c6<<"\tc7="<<c7<<"\tc8="<<c8<<"\tc9="<<c9<<std::endl;
//     out<<"Du="<<Du<<"\tDv="<<Dv<<"\tDw="<<Dw<<std::endl;
//     out<<"Fmax="<<Fmax<<"\tGmax="<<Gmax<<"\tHmax="<<Hmax<<std::endl;
//     out<<"P="<<P<<std::endl;
//     out<<"tmax="<<tmax<<std::endl;
//     out<<"Nx="<<Nx<<"\tNy="<<Ny<<std::endl;
//     out<<"dt="<<dt<<std::endl;
//     out<<"epsilon="<<epsilon<<std::endl;
//     out<<"S="<<S<<std::endl;

//     out<<"delta_obs="<<delta_obs<<std::endl;
//     out<<"result_folder="<<result_folder<<std::endl;
// 	out<<"********************************************"<<std::endl<<std::endl;
// 	return out;
// };

// void write_params(const Parameters &a, const std::string &name, const ios_base::openmode mode=ios_base::trunc)
// {
// 	ofstream file;
// 	file.open(name.c_str(),mode);
// 	file<<a;
// 	file.close();
// };

#endif
