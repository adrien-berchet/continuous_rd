#define BOOST_LOG_DYN_LINK 1

#include <boost/filesystem.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

#include "logger.hpp"

#include "rd_dynamics.hpp"
#include "class_parameters.hpp"
#include "definitions.hpp"
#include "hex_2d_lattice.hpp"
#include "initialization.hpp"

namespace ode = boost::numeric::odeint;

host_state_type simulate_rd(Parameters &params)
{
	// Get values from parameters
	const double &cu=params.cu, &cv=params.cv, &cw=params.cw;
	const double &c1=params.c1, &c2=params.c2, &c3=params.c3, &c4=params.c4, &c5=params.c5, &c6=params.c6, &c7=params.c7, &c8=params.c8, &c9=params.c9;
	const double &Du=params.Du, &Dv=params.Dv, &Dw=params.Dw;
	const double &Fmax=params.Fmax, &Gmax=params.Gmax, &Hmax=params.Hmax;

	const size_t &Nx = params.Nx, &Ny = params.Ny;
	const size_t N = Nx * Ny;
	const double &epsilon = params.epsilon;
	const value_type &dt = params.dt;

	// Create vectors of data: all variables are concatenated into one vector for simplicity
	// Create initial conditions and initial values on host
	host_state_type x_host( 3 * N, 0 );
	host_state_type Pxx_top( N, 1.0 ), Pxx_bot( N, 1.0 ), Pxx_left( N, 1.0 ), Pxx_right( N, 1.0 );
	if (params.gauss_std > 0)
	{
		gauss_init(x_host, params, Pxx_top, Pxx_bot, Pxx_left, Pxx_right);
	}
	else
	{
		random_hex_init(x_host, params, Pxx_top, Pxx_bot, Pxx_left, Pxx_right);
	}

	// Copy to device
	state_type x = x_host;

	// Create stepper
	ode::runge_kutta4< state_type , value_type , state_type , value_type > stepper;

	// Create the dynamic system function
	rd_dynamics sys(
		Nx, Ny, epsilon,
		cu, cv, cw,
		c1, c2, c3, c4, c5, c6, c7, c8, c9,
		Du, Dv, Dw,
		Fmax, Gmax, Hmax,
		Pxx_top, Pxx_bot, Pxx_left, Pxx_right
	);

	// Export neighbors if requested
	if(params.export_neighbors)
	{
		export_neighbors(sys, params);
	}

	// Create observer
	observer obs(params, N, sys.Pxx_top, sys.Pxx_bot, sys.Pxx_left, sys.Pxx_right);

	// Integrate
	// TODO: Add stoping criteria but Boost::ODEINT does not provide an easy way to do this. I think this should be done inside the observer to interrupt the integration when the criteria is satisfied. Another solution can be to just use do_step() manually.
	integrate_const( stepper , sys , x , 0.0 , params.tmax , dt , boost::ref(obs));

	// Export final state
	if (obs.last_t != params.tmax)
	{
		obs.last_t = -1;
		obs(x, params.tmax);
	}

	// Export results
	thrust::copy( x.begin() , x.end() , x_host.begin() );
	return x_host;
}

int main( int argc , char* argv[] )
{
	// Define and read the parameters
	Parameters params;
	params.read(argc, argv);
	std::cout<<params<<std::endl;

	// Initialize the logger
	generic::init_logger(params.log_level);

	// Create folders in which the results will be stored
	if(!boost::filesystem::is_directory(params.result_folder))
	{
		boost::filesystem::create_directories(params.result_folder);
	}
	if(!boost::filesystem::is_directory(params.result_folder + "/results"))
	{
		boost::filesystem::create_directories(params.result_folder + "/results");
	}

	// Export parameters used for the current simulation
	BOOST_LOG_TRIVIAL(info) << "Export parameters";
	params.write_parameters(params.result_folder + "/parameters_used.prm");

	// Run the simulation
	host_state_type result = simulate_rd(params);
}
