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
#include "dat_writer.hpp"

namespace ode = boost::numeric::odeint;


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
		const size_t &Nx_in, const size_t &Ny_in,
		const double &cu_in, const double &cv_in, const double &cw_in,
		const double &c1_in, const double &c2_in, const double &c3_in, const double &c4_in, const double &c5_in, const double &c6_in, const double &c7_in, const double &c8_in, const double &c9_in,
		const double &Du_in, const double &Dv_in, const double &Dw_in,
		const double &Fmax_in, const double &Gmax_in, const double &Hmax_in
	):
		N ( Nx_in * Ny_in ), Nx( Nx_in ), Ny ( Ny_in ),
		cu(cu_in), cv(cv_in), cw(cw_in),
		c1(c1_in), c2(c2_in), c3(c3_in), c4(c4_in), c5(c5_in), c6(c6_in), c7(c7_in), c8(c8_in), c9(c9_in),
		Du(Du_in), Dv(Dv_in), Dw(Dw_in),
		Fmax(Fmax_in), Gmax(Gmax_in), Hmax(Hmax_in),
		top( 3 * N ), bot( 3 * N ), left( 3 * N ), right( 3 * N )
	{
		// Define neighbors
		thrust::counting_iterator<size_t> counter( 0 );

		// Top neighbors
		thrust::copy( counter , counter + (N - Nx) , top.begin() + Nx ); // u component
		thrust::copy( counter + N , counter + (2 * N - Nx) , top.begin() + N + Nx ); // v component
		thrust::copy( counter + 2 * N , counter + (3 * N - Nx) , top.begin() + 2 * N + Nx ); // w component
		thrust::copy( counter + (N - Nx), counter + N , top.begin() ); // u component
		thrust::copy( counter + (2 * N - Nx), counter + 2 * N , top.begin() + N ); // v component
		thrust::copy( counter + (3 * N - Nx), counter + 3 * N , top.begin() + 2 * N); // w component

		// Bottom neighbors
		thrust::copy( counter + Nx , counter + N , bot.begin() ); // u component
		thrust::copy( counter + N + Nx , counter + 2 * N , bot.begin() + N ); // v component
		thrust::copy( counter + 2 * N + Nx , counter + 3 * N , bot.begin() + 2 * N ); // w component
		thrust::copy( counter, counter + Nx , bot.begin() + N - Nx ); // u component
		thrust::copy( counter + N, counter + N + Nx , bot.begin() + 2 * N - Nx ); // V component
		thrust::copy( counter + 2 * N, counter + 2 * N + Nx , bot.begin() + 3 * N - Nx ); // w component

		// Left neighbors
		thrust::copy( counter , counter + 3 * N - 1 , left.begin() + 1 );

		// Right neighbors
		thrust::copy( counter + 1 , counter + 3 * N , right.begin() );

		// Adjust left and right neighbors on sides
		for (int i = 0; i < Ny; ++i)
		{
			left[i * Nx] = (i + 1) * Nx - 1; // u component
			right[(i + 1) * Nx - 1] = i * Nx; // u component
			left[N + i * Nx] = N + (i + 1) * Nx - 1; // v component
			right[N + (i + 1) * Nx - 1] = N + i * Nx; // v component
			left[2 * N + i * Nx] = 2 * N + (i + 1) * Nx - 1; // w component
			right[2 * N + (i + 1) * Nx - 1] = 2 * N + i * Nx; // w component
		}
	}

	void operator() ( const state_type &x , state_type &dxdt , const value_type )
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
							thrust::make_permutation_iterator( x.begin(), top.begin() ),
							thrust::make_permutation_iterator( x.begin(), top.begin() + N ),
							thrust::make_permutation_iterator( x.begin(), top.begin() + 2 * N )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin(), bot.begin() ),
							thrust::make_permutation_iterator( x.begin(), bot.begin() + N ),
							thrust::make_permutation_iterator( x.begin(), bot.begin() + 2 * N )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin(), left.begin() ),
							thrust::make_permutation_iterator( x.begin(), left.begin() + N ),
							thrust::make_permutation_iterator( x.begin(), left.begin() + 2 * N )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin(), right.begin() ),
							thrust::make_permutation_iterator( x.begin(), right.begin() + N ),
							thrust::make_permutation_iterator( x.begin(), right.begin() + 2 * N )
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
							thrust::make_permutation_iterator( x.begin(), top.begin() + N ),
							thrust::make_permutation_iterator( x.begin(), top.begin() + 2 * N ),
							thrust::make_permutation_iterator( x.begin(), top.end() )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin(), bot.begin() + N ),
							thrust::make_permutation_iterator( x.begin(), bot.begin() + 2 * N ),
							thrust::make_permutation_iterator( x.begin(), bot.end() )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin(), left.begin() + N ),
							thrust::make_permutation_iterator( x.begin(), left.begin() + 2 * N ),
							thrust::make_permutation_iterator( x.begin(), left.end() )
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(
							thrust::make_permutation_iterator( x.begin(), right.begin() + N ),
							thrust::make_permutation_iterator( x.begin(), right.begin() + 2 * N ),
							thrust::make_permutation_iterator( x.begin(), right.end() )
					) ) ,
					thrust::make_constant_iterator( thrust::make_tuple(cu, cv, cw, Du, Dv, Dw) ),
					thrust::make_constant_iterator( thrust::make_tuple(c1, c2, c3, c4, c5, c6, c7, c8, c9) ),
					thrust::make_constant_iterator( thrust::make_tuple(Fmax, Gmax, Hmax) ),
					thrust::make_zip_iterator(
						thrust::make_tuple(dxdt.begin() + N, dxdt.begin() + 2 * N, dxdt.begin() + 3 * N)
					)
			) ),
			sys_functor()
		);
	}

	const index_vector_type& get_top() const {return this->top;}
	const index_vector_type& get_bot() const {return this->bot;}
	const index_vector_type& get_left() const {return this->left;}
	const index_vector_type& get_right() const {return this->right;}
	const size_t& get_N() const {return this->N;}

private:

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
    	// TODO: use params.delta_obs to skip some exports if they are too close from each other

    	// Format file name (zero padding to ensure that the file are always correctly sorted)
		std::ostringstream filename;
		filename << std::fixed << std::setprecision(precision) << std::setw(filename_length) << std::setfill('0') << t;

		// Create file
        generic::DatWriter data_file(params.result_folder + "/results/" + filename.str() + ".dat");

		// Write header
        data_file.write_header(std::to_string(t), params.Nx, params.Ny, "x", "y", "u", "v", "w", "P_sin_theta");

        // Write data
        // TODO: This is the slowest part of the code, try to improve it either by using binary files or by improving the iterators (but I/O are probably the main limitation)
        int num=0;
        for(
        	auto i=thrust::make_zip_iterator(thrust::make_tuple(
				state.begin(),
				state.begin() + N,
				state.begin() + 2 * N,
				state.begin() + 3 * N
			) );
			i != thrust::make_zip_iterator(thrust::make_tuple(
				state.begin() + N,
				state.begin() + 2 * N,
				state.begin() + 3 * N,
				state.begin() + 4 * N
			) );
			++i
		)
        {
        	const double x = (num % params.Nx) * params.epsilon;
        	const double y = int(num / params.Nx) * params.epsilon;
        	data_file.write_row(x, y, thrust::get<0>(*i), thrust::get<1>(*i), thrust::get<2>(*i), thrust::get<3>(*i));
        	++num;
        }
    }
};

/**
 * \brief Random initialization
*/
void random_init(std::vector< value_type > &state, const Parameters &params)
{}

/**
 * \brief Gaussian initialization. Only used for validation.
*/
void gauss_init(std::vector< value_type > &state, const Parameters &params)
{
	const size_t &Nx = params.Nx, &Ny = params.Ny;
	const size_t N = Nx * Ny;

	const double &epsilon = params.epsilon;
	const double &sigma = params.gauss_std;

	const double x_center = (Nx / 2.0) * epsilon;
	const double y_center = (Ny / 2.0) * epsilon;

	std::cout<<"Gaussian initialization:"<<std::endl;
	std::cout<<"\t x_center = "<<x_center<<std::endl;
	std::cout<<"\t y_center = "<<y_center<<std::endl;
	std::cout<<"\t epsilon = "<<epsilon<<std::endl;
	std::cout<<"\t std = "<<sigma<<std::endl;

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
    	const double x = (num % Nx) * epsilon;
    	const double y = int(num / Nx) * epsilon;
		const double r = std::sqrt(std::pow(x - x_center, 2.) + std::pow(y - y_center, 2.));
		const double C = std::exp(-std::pow(r, 2.) / (2. * std::pow(sigma, 2.)));
    	thrust::get<0>(*i) = C;
    	thrust::get<1>(*i) = C;
    	thrust::get<2>(*i) = C;
    	++num;
    }
}

void export_neighbors(const rd_dynamics &sys, const Parameters &params)
{
	const size_t &N=sys.get_N();

	// Create file
    generic::DatWriter data_file(params.result_folder + "/neighbors.dat");

	// Write header
    data_file.write_header("Neighbors", params.Nx, params.Ny, "x", "y", "num", "top_u", "top_v", "top_w", "bot_u", "bot_v", "bot_w", "left_u", "left_v", "left_w", "right_u", "right_v", "right_w");

    // Write data
    int num=0;
    for(
    	auto i=thrust::make_zip_iterator(
			thrust::make_tuple(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						sys.get_top().begin(),
						sys.get_top().begin() + N,
						sys.get_top().begin() + 2 * N
				) ) ,
				thrust::make_zip_iterator(
					thrust::make_tuple(
						sys.get_bot().begin() ,
						sys.get_bot().begin() + N ,
						sys.get_bot().begin() + 2 * N
				) ) ,
				thrust::make_zip_iterator(
					thrust::make_tuple(
						sys.get_left().begin() ,
						sys.get_left().begin() + N ,
						sys.get_left().begin() + 2 * N
				) ) ,
				thrust::make_zip_iterator(
					thrust::make_tuple(
						sys.get_right().begin() ,
						sys.get_right().begin() + N ,
						sys.get_right().begin() + 2 * N
				) )
		) );
		i != thrust::make_zip_iterator(
			thrust::make_tuple(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						sys.get_top().begin() + N ,
						sys.get_top().begin() + 2 * N ,
						sys.get_top().end()
				) ) ,
				thrust::make_zip_iterator(
					thrust::make_tuple(
						sys.get_bot().begin() + N ,
						sys.get_bot().begin() + 2 * N ,
						sys.get_bot().end()
				) ) ,
				thrust::make_zip_iterator(
					thrust::make_tuple(
						sys.get_left().begin() + N ,
						sys.get_left().begin() + 2 * N ,
						sys.get_left().end()
				) ) ,
				thrust::make_zip_iterator(
					thrust::make_tuple(
						sys.get_right().begin() + N ,
						sys.get_right().begin() + 2 * N ,
						sys.get_right().end()
				) )
		) );
		++i
	)
    {
    	const double x = (num % params.Nx) * params.epsilon;
    	const double y = int(num / params.Nx) * params.epsilon;
    	data_file.write_row(
    		x, y, num,
    		thrust::get<0>(thrust::get<0>(*i)), thrust::get<1>(thrust::get<0>(*i)), thrust::get<2>(thrust::get<0>(*i)),
    		thrust::get<0>(thrust::get<1>(*i)), thrust::get<1>(thrust::get<1>(*i)), thrust::get<2>(thrust::get<1>(*i)),
    		thrust::get<0>(thrust::get<2>(*i)), thrust::get<1>(thrust::get<2>(*i)), thrust::get<2>(thrust::get<2>(*i)),
    		thrust::get<0>(thrust::get<3>(*i)), thrust::get<1>(thrust::get<3>(*i)), thrust::get<2>(thrust::get<3>(*i))
		);
    	++num;
    }
}

std::vector<value_type> simulate_rd(Parameters &params)
{
	// Get values from parameters
	const double &cu=params.cu, &cv=params.cv, &cw=params.cw;
	const double &c1=params.c1, &c2=params.c2, &c3=params.c3, &c4=params.c4, &c5=params.c5, &c6=params.c6, &c7=params.c7, &c8=params.c8, &c9=params.c9;
	const double &Du=params.Du, &Dv=params.Dv, &Dw=params.Dw;
	const double &Fmax=params.Fmax, &Gmax=params.Gmax, &Hmax=params.Hmax;

	const size_t &Nx = params.Nx, &Ny = params.Ny;
	const size_t N = Nx * Ny;
	const value_type &dt = params.dt;

	// Create vectors of data: all variables are concatenated into one vector for simplicity
	// Create initial conditions and initial values on host
	std::vector< value_type > x_host( 4 * N, 0 );
    std::fill(x_host.begin() + 3 * N, x_host.end(), 1.0); // Set P sin(theta) = 1 everywhere
	if (params.gauss_std > 0)
	{
		gauss_init(x_host, params);
	}
	else
	{
		random_init(x_host, params);
	}

	// Copy to device
	state_type x = x_host;

	// Create stepper
	ode::runge_kutta4< state_type , value_type , state_type , value_type > stepper;

	// Create phase oscillator system function
	rd_dynamics sys(
		Nx, Ny,
		cu, cv, cw,
		c1, c2, c3, c4, c5, c6, c7, c8, c9,
		Du, Dv, Dw,
		Fmax, Gmax, Hmax
	);

	// Export neighbors
	if(params.export_neighbors)
	{
		export_neighbors(sys, params);
	}

	// Create observer
	observer obs(params, N);

	// Integrate
	// TODO: Add stoping criteria but Boost::ODEINT does not provide an easy way to do this. I think this should be done inside the observer to interrupt the integration when the criteria is satisfied. Another solution can be to just use do_step() manually.
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
	std::cout<<params<<std::endl;

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
	params.write_parameters(params.result_folder + "/parameters_used.prm");

	// Run the simulation
	std::vector<value_type> result = simulate_rd(params);
}
