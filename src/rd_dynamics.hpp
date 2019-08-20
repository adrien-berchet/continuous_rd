/**
 * \file rd_dynamics.hpp
 *
 * \brief Provides a class and some tools to compute the RD dnyamics using
 * Boost::ODEINT and CUDA.
 *
*/
#ifndef RD_DYNAMICS_HPP
#define RD_DYNAMICS_HPP

#include <iomanip>

#include <thrust/device_vector.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#include "dat_writer.hpp"
#include "logger.hpp"

#include "class_parameters.hpp"
#include "definitions.hpp"


/**
 * \brief Define the dynamic system
*/
class rd_dynamics
{

public:

	/**
	 * \brief Structure that define the functor that compute dx/dt at a given time step.
	*/
	struct sys_functor
	{
		/**
		 * \brief Functor that compute dx/dt at a given time step.
		*/
		template< class Tuple >
		__host__ __device__
		void operator()( Tuple t )  // This functor works on tuples of values
		{
			// Get current values
			const value_type u = thrust::get<0>(thrust::get<0>(t));
			const value_type v = thrust::get<1>(thrust::get<0>(t));
			const value_type w = thrust::get<2>(thrust::get<0>(t));
			const double epsilon = thrust::get<3>(thrust::get<0>(t));

			// Get P sin(theta) for each direction
			const value_type P_sin_theta_top = thrust::get<0>(thrust::get<8>(t));
			const value_type P_sin_theta_bot = thrust::get<1>(thrust::get<8>(t));
			const value_type P_sin_theta_left = thrust::get<2>(thrust::get<8>(t));
			const value_type P_sin_theta_right = thrust::get<3>(thrust::get<8>(t));

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

			const value_type lapl_u = (
				  P_sin_theta_top * (u_top - u)
				+ P_sin_theta_bot * (u_bot - u)
				+ P_sin_theta_left * (u_left - u)
				+ P_sin_theta_right * (u_right - u));
			const value_type lapl_v = (
				  P_sin_theta_top * (v_top - v)
				+ P_sin_theta_bot * (v_bot - v)
				+ P_sin_theta_left * (v_left - v)
				+ P_sin_theta_right * (v_right - v));
			const value_type lapl_w = (
				  P_sin_theta_top * (w_top - w)
				+ P_sin_theta_bot * (w_bot - w)
				+ P_sin_theta_left * (w_left - w)
				+ P_sin_theta_right * (w_right - w));

			// The dynamical equation
			thrust::get<0>(thrust::get<9>(t)) = F - cu * u + Du * lapl_u;
			thrust::get<1>(thrust::get<9>(t)) = G - cv * v + Dv * lapl_v;
			thrust::get<2>(thrust::get<9>(t)) = H - cw * w + Dw * lapl_w;
		}
	};

	/**
	 * \brief Constructor. Initialize all attributes and the neighbors used to cumpute the Laplacian.
	*/
	rd_dynamics(
		const size_t &Nx_in, const size_t &Ny_in, const double &epsilon_in,
		const double &cu_in, const double &cv_in, const double &cw_in,
		const double &c1_in, const double &c2_in, const double &c3_in, const double &c4_in, const double &c5_in, const double &c6_in, const double &c7_in, const double &c8_in, const double &c9_in,
		const double &Du_in, const double &Dv_in, const double &Dw_in,
		const double &Fmax_in, const double &Gmax_in, const double &Hmax_in,
		host_state_type &Pxx_top_in, host_state_type &Pxx_bot_in, host_state_type &Pxx_left_in, host_state_type &Pxx_right_in
	):
		N ( Nx_in * Ny_in ), Nx( Nx_in ), Ny ( Ny_in ), epsilon( epsilon_in ),
		cu(cu_in), cv(cv_in), cw(cw_in),
		c1(c1_in), c2(c2_in), c3(c3_in), c4(c4_in), c5(c5_in), c6(c6_in), c7(c7_in), c8(c8_in), c9(c9_in),
		Du(Du_in), Dv(Dv_in), Dw(Dw_in),
		Fmax(Fmax_in), Gmax(Gmax_in), Hmax(Hmax_in),
		top( 3 * N ), bot( 3 * N ), left( 3 * N ), right( 3 * N ),
		Pxx_top( Pxx_top_in ), Pxx_bot( Pxx_bot_in ), Pxx_left( Pxx_left_in ), Pxx_right( Pxx_right_in )
	{
		BOOST_LOG_TRIVIAL(debug) << "RD dynamics initialization";

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

	/**
	 * \brief Functor called by Boost::ODEINT at each time step.
	*/
	void operator() ( const state_type &x , state_type &dxdt , const value_type &t)
	{
		BOOST_LOG_TRIVIAL(debug) << "Compute dynamics for t="<<t;
		thrust::for_each(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					thrust::make_zip_iterator(
						thrust::make_tuple(
							x.begin() ,
							x.begin() + N,
							x.begin() + 2 * N,
							thrust::make_constant_iterator(epsilon)
					) ) ,
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
						thrust::make_tuple(
							Pxx_top.begin(),
							Pxx_bot.begin(),
							Pxx_left.begin(),
							Pxx_right.begin()
					) ) ,
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
							thrust::make_constant_iterator(epsilon)
					) ) ,
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
						thrust::make_tuple(
							Pxx_top.end(),
							Pxx_bot.end(),
							Pxx_left.end(),
							Pxx_right.end()
					) ) ,
					thrust::make_zip_iterator(
						thrust::make_tuple(dxdt.begin() + N, dxdt.begin() + 2 * N, dxdt.begin() + 3 * N)
					)
			) ),
			sys_functor()
		);
	}

	/** \brief Return top neighbors. */
	const index_vector_type& get_top() const {return this->top;}
	/** \brief Return bot neighbors. */
	const index_vector_type& get_bot() const {return this->bot;}
	/** \brief Return left neighbors. */
	const index_vector_type& get_left() const {return this->left;}
	/** \brief Return right neighbors. */
	const index_vector_type& get_right() const {return this->right;}
	/** \brief Return the number of elements in the regular grid. */
	const size_t& get_N() const {return this->N;}

private:
	/**@{ */
	const size_t N, Nx, Ny;
	const double cu, cv, cw;
	const double c1, c2, c3, c4, c5, c6, c7, c8, c9;
	const double Du, Dv, Dw;
	const double Fmax, Gmax, Hmax;
	const double epsilon;
	index_vector_type top, bot, left, right;
	/**@}*/

public:
	/**@{ */
	const state_type Pxx_top, Pxx_bot, Pxx_left, Pxx_right;
	/**@}*/
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
	/**
	 * \cond DOXYGEN_IGNORE
	*/
	const Parameters &params;
	const size_t N;
	const size_t filename_length;
	const size_t precision;
	const state_type &Pxx_top, &Pxx_bot, &Pxx_left, &Pxx_right;
	double last_t;
	/**
	 * \endcond
	*/


	/**
	 * \brief Constructor.
	*/
	observer(
		const Parameters &params_in, const size_t &N_in,
		const state_type &Pxx_top_in,
		const state_type &Pxx_bot_in,
		const state_type &Pxx_left_in,
		const state_type &Pxx_right_in
	):
		params( params_in ), N( N_in ),
		filename_length( number_length(params.tmax, params.dt) ),
		precision( number_precision(params.dt) ),
		Pxx_top( Pxx_top_in ), Pxx_bot( Pxx_bot_in ), Pxx_left( Pxx_left_in ), Pxx_right( Pxx_right_in )
	{
		BOOST_LOG_TRIVIAL(debug) << "Observer initialization";
		last_t = - 1.01 * params_in.delta_obs;
	}

	/**
	 * \brief Functor called by Boost::ODEINT at each time step.
	*/
	template< class State >
	void operator()( const State &state , value_type t )
	{
		// Skip some exports if they are too close from each other
		if(t - last_t < params.delta_obs)
		{
			BOOST_LOG_TRIVIAL(debug) << "Skip export results for t="<<t;
			return;
		}
		else
		{
			if (last_t < 0)
			{
				last_t = t;
			}
			else
			{
				last_t += params.delta_obs;
			}
		}

		// Format file name (zero padding to ensure that the file are always correctly sorted)
		std::ostringstream filename;
		filename << std::fixed << std::setprecision(precision) << std::setw(filename_length) << std::setfill('0') << params.result_folder << "/results/" << t << ".dat";

		BOOST_LOG_TRIVIAL(info) << "Export results for t="<<t<<" in "<<filename.str();

		// Create file
		generic::DatWriter data_file(filename.str());

		// Write header
		data_file.write_header(std::to_string(t), params.Nx, params.Ny, "x", "y", "u", "v", "w", "P_sin_theta_top", "P_sin_theta_bot", "P_sin_theta_left", "P_sin_theta_right");

		// Transfert data to host
		host_state_type tmp = state;

		// Write data
		// TODO: This is the slowest part of the code, try to improve it either by using binary files or by improving the iterators (but I/O are probably the main limitation)
		int num=0;
		for(
			auto i=thrust::make_zip_iterator(thrust::make_tuple(
				tmp.begin(),
				tmp.begin() + N,
				tmp.begin() + 2 * N,
				Pxx_top.begin(),
				Pxx_bot.begin(),
				Pxx_left.begin(),
				Pxx_right.begin()
			) );
			i != thrust::make_zip_iterator(thrust::make_tuple(
				tmp.begin() + N,
				tmp.begin() + 2 * N,
				tmp.begin() + 3 * N,
				Pxx_top.end(),
				Pxx_bot.end(),
				Pxx_left.end(),
				Pxx_right.end()
			) );
			++i
		)
		{
			data_file.write_row(
				get_x(num, params.Nx, params.epsilon),
				get_y(num, params.Nx, params.epsilon),
				thrust::get<0>(*i),
				thrust::get<1>(*i),
				thrust::get<2>(*i),
				thrust::get<3>(*i),
				thrust::get<4>(*i),
				thrust::get<5>(*i),
				thrust::get<6>(*i)
			);
			++num;
		}
	}
};

/**
 * \brief Export neighbors indices used to compute the laplacian. Only used for validation.
*/
void export_neighbors(const rd_dynamics &sys, const Parameters &params)
{
	BOOST_LOG_TRIVIAL(info) << "Export neighbors";
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
		const double x = get_x(num, params.Nx, params.epsilon);
		const double y = get_y(num, params.Nx, params.epsilon);
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

#endif // RD_DYNAMICS_HPP
