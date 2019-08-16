#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <thrust/device_vector.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#include <boost/filesystem.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

#include "class_parameters.hpp"
#include "csv_writer.hpp"

using namespace std;

using namespace boost::numeric::odeint;


// Change this to float if your device does not support double computation
typedef double value_type;


#ifdef WITH_GPU
typedef thrust::device_vector< value_type > state_type;
typedef thrust::device_vector< size_t > index_vector_type;
#else
typedef thrust::host_vector< value_type > state_type;
typedef thrust::host_vector< size_t > index_vector_type;
#endif


/**
 * \brief Define the dynamic system
*/
class rd_dynamics
{

public:

	struct sys_functor
	{
		template< class Tuple >
		__host__ __device__
		void operator()( Tuple t )  // This functor works on tuples of values
		{
			// Get current values
			const value_type u = thrust::get<0>(thrust::get<0>(t));
			const value_type v = thrust::get<1>(thrust::get<0>(t));
			const value_type w = thrust::get<2>(thrust::get<0>(t));

			// Get P sin(theta)
			const value_type P_sin_theta = thrust::get<3>(thrust::get<0>(t));

			// Get neighbors
			const value_type u_top = thrust::get<0>(thrust::get<1>(t));  // Top neighbor for u
			const value_type v_top = thrust::get<1>(thrust::get<1>(t));  // Top neighbor for v
			const value_type w_top = thrust::get<2>(thrust::get<1>(t));  // Top neighbor for w
			const value_type u_bot = thrust::get<0>(thrust::get<2>(t));  // Bottom neighbor for u
			const value_type v_bot = thrust::get<1>(thrust::get<2>(t));  // Bottom neighbor for v
			const value_type w_bot = thrust::get<2>(thrust::get<2>(t));  // Bottom neighbor for w
			const value_type u_left = thrust::get<0>(thrust::get<3>(t));  // Left neighbor for u
			const value_type v_left = thrust::get<1>(thrust::get<3>(t));  // Left neighbor for v
			const value_type w_left = thrust::get<2>(thrust::get<3>(t));  // Left neighbor for w
			const value_type u_right = thrust::get<0>(thrust::get<4>(t)); // Right neighbor for u
			const value_type v_right = thrust::get<1>(thrust::get<4>(t)); // Right neighbor for v
			const value_type w_right = thrust::get<2>(thrust::get<4>(t)); // Right neighbor for w

			// Get cu, cv, cw
			const value_type cu = thrust::get<0>(thrust::get<5>(t));
			const value_type cv = thrust::get<1>(thrust::get<5>(t));
			const value_type cw = thrust::get<2>(thrust::get<5>(t));

			// Get Du, Dv, Dw
			const value_type Du = thrust::get<3>(thrust::get<5>(t));
			const value_type Dv = thrust::get<4>(thrust::get<5>(t));
			const value_type Dw = thrust::get<5>(thrust::get<5>(t));

			// Get c1, c2, ..., c9
			const value_type c1 = thrust::get<0>(thrust::get<6>(t));
			const value_type c2 = thrust::get<1>(thrust::get<6>(t));
			const value_type c3 = thrust::get<2>(thrust::get<6>(t));
			const value_type c4 = thrust::get<3>(thrust::get<6>(t));
			const value_type c5 = thrust::get<4>(thrust::get<6>(t));
			const value_type c6 = thrust::get<5>(thrust::get<6>(t));
			const value_type c7 = thrust::get<6>(thrust::get<6>(t));
			const value_type c8 = thrust::get<7>(thrust::get<6>(t));
			const value_type c9 = thrust::get<8>(thrust::get<6>(t));

			// Get Fmax, Gmax, Hmax
			const value_type Fmax = thrust::get<0>(thrust::get<7>(t));
			const value_type Gmax = thrust::get<1>(thrust::get<7>(t));
			const value_type Hmax = thrust::get<2>(thrust::get<7>(t));

			// Compute each term for each component
			const value_type F_cond = c1 * v + c2 * w + c3;
			const value_type G_cond = c4 * u + c5 * w + c6;
			const value_type H_cond = c7 * u + c8 * v + c9;

			const value_type F = thrust::max(0.0, thrust::min(Fmax, F_cond));
			const value_type G = thrust::max(0.0, thrust::min(Gmax, G_cond));
			const value_type H = thrust::max(0.0, thrust::min(Hmax, H_cond));

			const value_type lapl_u = u_top + u_bot + u_left + u_right - 4 * u;
			const value_type lapl_v = v_top + v_bot + v_left + v_right - 4 * v;
			const value_type lapl_w = w_top + w_bot + w_left + w_right - 4 * w;

			// The dynamical equation
			thrust::get<0>(thrust::get<8>(t)) = F - cu * u + Du * lapl_u * P_sin_theta;
			thrust::get<1>(thrust::get<8>(t)) = G - cv * v + Dv * lapl_v * P_sin_theta;
			thrust::get<2>(thrust::get<8>(t)) = H - cw * w + Dw * lapl_w * P_sin_theta;
		}
	};

	rd_dynamics(
		const state_type &init,
		const size_t &Nx_in, const size_t &Ny_in,
		const double &cu_in, const double &cv_in, const double &cw_in,
		const double &c1_in, const double &c2_in, const double &c3_in, const double &c4_in, const double &c5_in, const double &c6_in, const double &c7_in, const double &c8_in, const double &c9_in,
		const double &Du_in, const double &Dv_in, const double &Dw_in,
		const double &Fmax_in, const double &Gmax_in, const double &Hmax_in
	):
		state_n (init) , N ( Nx_in * Ny_in ), Nx( Nx_in ), Ny ( Ny_in ),
		cu(cu_in), cv(cv_in), cw(cw_in),
		c1(c1_in), c2(c2_in), c3(c3_in), c4(c4_in), c5(c5_in), c6(c6_in), c7(c7_in), c8(c8_in), c9(c9_in),
		Du(Du_in), Dv(Dv_in), Dw(Dw_in),
		Fmax(Fmax_in), Gmax(Gmax_in), Hmax(Hmax_in),
		top( init.size() ), bot( init.size() ), left( init.size() ), right( init.size() )
	{
		// Define neighbours
		thrust::counting_iterator<size_t> counter( 0 );

		// Top neighbours
		thrust::copy( counter , counter+(N-Nx) , top.begin()+Nx );
		thrust::copy( counter+(N-Nx), counter+N , top.begin() );

		// Bottom neighbours
		thrust::copy( counter+Nx , counter+N , bot.begin() );
		thrust::copy( counter, counter+Nx , bot.begin()+N-Nx );

		// Left neighbours
		thrust::copy( counter , counter+N-1 , left.begin()+1 );

		// Right neighbours
		thrust::copy( counter+1 , counter+N , right.begin() );

		// Adjust left and right neighbours on sides
		for (int i = 0; i < Ny; ++i)
		{
			left[i * Nx] = i * Nx - 1;
			right[(i + 1) * Nx - 1] = i * Nx;
		}
	}

	void operator() ( const state_type &x , state_type &dxdt , const value_type dt )
	{
		thrust::for_each(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					thrust::make_zip_iterator(
						thrust::make_tuple(
							x.begin() ,
							x.begin() + N,
							x.begin() + 2 * N,
							x.begin() + 3 * N
					) ),
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin() , top.begin() ) ,
							thrust::make_permutation_iterator( x.begin() + N , top.begin() + N ) ,
							thrust::make_permutation_iterator( x.begin() + 2 * N , top.begin() + 2 * N )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin() , bot.begin() ) ,
							thrust::make_permutation_iterator( x.begin() + N , bot.begin() + N ) ,
							thrust::make_permutation_iterator( x.begin() + 2 * N , bot.begin() + 2 * N )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin() , left.begin() ) ,
							thrust::make_permutation_iterator( x.begin() + N , left.begin() + N ) ,
							thrust::make_permutation_iterator( x.begin() + 2 * N , left.begin() + 2 * N )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin() , right.begin() ) ,
							thrust::make_permutation_iterator( x.begin() + N , right.begin() + N ) ,
							thrust::make_permutation_iterator( x.begin() + 2 * N , right.begin() + 2 * N )
					) ) ,
					thrust::make_constant_iterator( thrust::make_tuple(cu, cv, cw, Du, Dv, Dw) ),
					thrust::make_constant_iterator( thrust::make_tuple(c1, c2, c3, c4, c5, c6, c7, c8, c9) ),
					thrust::make_constant_iterator( thrust::make_tuple(Fmax, Gmax, Hmax) ),
					thrust::make_zip_iterator(
						thrust::make_tuple(dxdt.begin(), dxdt.begin() + N, dxdt.begin() + 2 * N)
					)
			) ),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					thrust::make_zip_iterator(
						thrust::make_tuple(
							x.begin() + N,
							x.begin() + 2 * N,
							x.begin() + 3 * N,
							x.end()
					) ),
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin() + N , top.begin() + N ) ,
							thrust::make_permutation_iterator( x.begin() + 2 * N , top.begin() + 2 * N ),
							thrust::make_permutation_iterator( x.end() , top.end() )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin() + N , bot.begin() + N ) ,
							thrust::make_permutation_iterator( x.begin() + 2 * N , bot.begin() + 2 * N ),
							thrust::make_permutation_iterator( x.end() , bot.end() )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin() + N , left.begin() + N ) ,
							thrust::make_permutation_iterator( x.begin() + 2 * N , left.begin() + 2 * N ),
							thrust::make_permutation_iterator( x.end() , left.end() )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin() + N , right.begin() + N ) ,
							thrust::make_permutation_iterator( x.begin() + 2 * N , right.begin() + 2 * N ),
							thrust::make_permutation_iterator( x.end() , right.end() )
					) ) ,
					thrust::make_constant_iterator( thrust::make_tuple(cu, cv, cw, Du, Dv, Dw) ),
					thrust::make_constant_iterator( thrust::make_tuple(c1, c2, c3, c4, c5, c6, c7, c8, c9) ),
					thrust::make_constant_iterator( thrust::make_tuple(Fmax, Gmax, Hmax) ),
					thrust::make_zip_iterator(
						thrust::make_tuple(dxdt.begin() + N, dxdt.begin() + 2 * N, dxdt.end())
					)
			) ),
			sys_functor()
		);
	}

private:

	const state_type &state_n;
	const size_t N, Nx, Ny;
	const double cu, cv, cw;
	const double c1, c2, c3, c4, c5, c6, c7, c8, c9;
	const double Du, Dv, Dw;
	const double Fmax, Gmax, Hmax;
	index_vector_type top, bot, left, right;
};

/**
 * \brief Compute the max length of the file names
*/
template<typename T>
size_t number_length(T &tmax, T &dt)
{
	std::ostringstream tmp;
	double int_part;
	double decimal_part = std::modf(tmax, &int_part);
	tmp << int_part + dt;
	return tmp.str().size();
}

/**
 * \brief Compute the max length of the decimal part of the file names
*/
template<typename T>
size_t number_precision(T &dt)
{
	std::ostringstream tmp;
	tmp << dt;
	return std::max(size_t(3), tmp.str().size()) - 2;
}

/**
 * \brief Define the observer used to export the results
*/
struct observer
{
    const Parameters &params;
    const size_t N;
    const size_t filename_length;
    const size_t precision;

    observer( const Parameters &params_in, const size_t &N_in ) : params( params_in ), N( N_in ), filename_length( number_length(params.tmax, params.dt) ), precision( number_precision(params.dt) ) {}

    template< class State >
    void operator()( const State &state , value_type t )
    {
    	// Format file name (zero padding to ensure that the file are always correctly sorted)
		std::ostringstream filename;
		filename << std::fixed << std::setprecision(precision) << std::setw(filename_length) << std::setfill('0') << t;

		// Create file
        generic::CsvWriter csv_file(params.result_folder + "/results/" + filename.str() + ".dat");

		// Write header
        csv_file.write_header("x", "y", "u", "v", "w");

        // Write data
        // TODO: This is the slowest part of the code, try to improve it (I/O are probably limitating but it may be possible to improve the iterators)
        int num=0;
        for(
        	auto i=thrust::make_zip_iterator(thrust::make_tuple(
				state.begin(),
				state.begin() + N,
				state.begin() + 2 * N
			) );
			i != thrust::make_zip_iterator(thrust::make_tuple(
				state.begin() + N,
				state.begin() + 2 * N,
				state.begin() + 3 * N
			) );
			++i
		)
        {
        	const double x = (num % params.Nx) * params.epsilon;
        	const double y = int(num / params.Nx) * params.epsilon;
        	csv_file.write_row(x, y, thrust::get<0>(*i), thrust::get<1>(*i), thrust::get<2>(*i));
        	++num;
        }
    }

};

std::vector<value_type> simulate_rd(Parameters &params)
{
	// Get values from parameters
	const double &cu=params.cu, &cv=params.cv, &cw=params.cw;
	const double &c1=params.c1, &c2=params.c2, &c3=params.c3, &c4=params.c4, &c5=params.c5, &c6=params.c6, &c7=params.c7, &c8=params.c8, &c9=params.c9;
	const double &Du=params.Du, &Dv=params.Dv, &Dw=params.Dw;
	const double &Fmax=params.Fmax, &Gmax=params.Gmax, &Hmax=params.Hmax;
	const double &P=params.P;

	const size_t &Nx = params.Nx, &Ny = params.Ny;
	const size_t N = Nx * Ny;
	const value_type &dt = params.dt;

	// Create vectors of data: all variables are concatenated into one vector for simplicity
	// Create initial conditions and initial values on host
	vector< value_type > x_host( 4 * N, 0 );
	vector< value_type > init_host( 4 * N, 0 );
	for( size_t i=0 ; i<(3 * N) ; ++i )
	{
		x_host[i] = 2.0 * drand48();
		init_host[i] = ( 4 * N - i ); // decreasing frequencies
	}
	for( size_t i=3 * N ; i<(4 * N) ; ++i )
	{
		x_host[i] = 1 + P;
		init_host[i] = 1 + P;
	}

	// Copy to device
	state_type x = x_host;
	state_type init = init_host;

	// Create stepper
	runge_kutta4< state_type , value_type , state_type , value_type > stepper;

	// Create phase oscillator system function
	rd_dynamics sys(
		init,
		Nx, Ny,
		cu, cv, cw,
		c1, c2, c3, c4, c5, c6, c7, c8, c9,
		Du, Dv, Dw,
		Fmax, Gmax, Hmax
	);

	// Create observer
	observer obs(params, N);

	// Integrate
	integrate_const( stepper , sys , x , 0.0 , params.tmax , dt , boost::ref(obs));

	// Export results
	thrust::copy( x.begin() , x.end() , x_host.begin() );
	return x_host;
}

int main( int argc , char* argv[] )
{
	// Define and read the parameters
	Parameters params;
	params.read(argc, argv);
	cout<<params<<endl;

	// Create folders in which the results will be stored
	if(!boost::filesystem::is_directory(params.result_folder))
	{
		boost::filesystem::create_directories(params.result_folder);
	}
	if(!boost::filesystem::is_directory(params.result_folder + "/results"))
	{
		boost::filesystem::create_directories(params.result_folder + "/results");
	}
	params.write_parameters(params.result_folder + "/parameters_used.prm");

	// Run the simulation
	vector<value_type> result = simulate_rd(params);
}
