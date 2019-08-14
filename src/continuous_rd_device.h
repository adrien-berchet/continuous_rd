#include <iostream>
#include <cmath>

#include <thrust/device_vector.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

// #include "class_parameters.hpp"

using namespace std;

using namespace boost::numeric::odeint;


// change this to float if your device does not support double computation
typedef double value_type;


#ifdef WITH_GPU
typedef thrust::device_vector< value_type > state_type;
typedef thrust::device_vector< size_t > index_vector_type;
#else
typedef thrust::host_vector< value_type > state_type;
typedef thrust::host_vector< size_t > index_vector_type;
#endif


class rd_dynamics
{

public:

	struct sys_functor
	{
		template< class Tuple >
		__host__ __device__
		void operator()( Tuple t )  // this functor works on tuples of values
		// {
		// 	// get current values
		// 	const value_type u = thrust::get<0>(thrust::get<0>(t));
		// 	const value_type v = thrust::get<1>(thrust::get<0>(t));
		// 	const value_type w = thrust::get<2>(thrust::get<0>(t));

		// 	// get P sin(theta)
		// 	const value_type P_sin_theta = thrust::get<3>(thrust::get<0>(t));

		// 	// get neighbors
		// 	const value_type u_top = thrust::get<0>(thrust::get<1>(t));  // top neighbor for u
		// 	const value_type v_top = thrust::get<1>(thrust::get<1>(t));  // top neighbor for v
		// 	const value_type w_top = thrust::get<2>(thrust::get<1>(t));  // top neighbor for w
		// 	const value_type u_bot = thrust::get<0>(thrust::get<2>(t));  // bottom neighbor for u
		// 	const value_type v_bot = thrust::get<1>(thrust::get<2>(t));  // bottom neighbor for v
		// 	const value_type w_bot = thrust::get<2>(thrust::get<2>(t));  // bottom neighbor for w
		// 	const value_type u_left = thrust::get<0>(thrust::get<3>(t));  // left neighbor for u
		// 	const value_type v_left = thrust::get<1>(thrust::get<3>(t));  // left neighbor for v
		// 	const value_type w_left = thrust::get<2>(thrust::get<3>(t));  // left neighbor for w
		// 	const value_type u_right = thrust::get<0>(thrust::get<4>(t)); // right neighbor for u
		// 	const value_type v_right = thrust::get<1>(thrust::get<4>(t)); // right neighbor for v
		// 	const value_type w_right = thrust::get<2>(thrust::get<4>(t)); // right neighbor for w

		// 	// get cu, cv, cw
		// 	const value_type cu = thrust::get<0>(thrust::get<5>(t));
		// 	const value_type cv = thrust::get<1>(thrust::get<5>(t));
		// 	const value_type cw = thrust::get<2>(thrust::get<5>(t));

		// 	// get Du, Dv, Dw
		// 	const value_type Du = thrust::get<3>(thrust::get<5>(t));
		// 	const value_type Dv = thrust::get<4>(thrust::get<5>(t));
		// 	const value_type Dw = thrust::get<5>(thrust::get<5>(t));

		// 	// get c1, c2, ..., c9
		// 	const value_type c1 = thrust::get<0>(thrust::get<6>(t));
		// 	const value_type c2 = thrust::get<1>(thrust::get<6>(t));
		// 	const value_type c3 = thrust::get<2>(thrust::get<6>(t));
		// 	const value_type c4 = thrust::get<3>(thrust::get<6>(t));
		// 	const value_type c5 = thrust::get<4>(thrust::get<6>(t));
		// 	const value_type c6 = thrust::get<5>(thrust::get<6>(t));
		// 	const value_type c7 = thrust::get<6>(thrust::get<6>(t));
		// 	const value_type c8 = thrust::get<7>(thrust::get<6>(t));
		// 	const value_type c9 = thrust::get<8>(thrust::get<6>(t));

		// 	// get Fmax, Gmax, Hmax
		// 	const value_type Fmax = thrust::get<0>(thrust::get<7>(t));
		// 	const value_type Gmax = thrust::get<1>(thrust::get<7>(t));
		// 	const value_type Hmax = thrust::get<2>(thrust::get<7>(t));

		// 	// compute each term for each component
		// 	const value_type F_cond = c1 * v + c2 * w + c3;
		// 	const value_type G_cond = c4 * u + c5 * w + c6;
		// 	const value_type H_cond = c7 * u + c8 * v + c9;

		// 	const value_type F = thrust::max(0.0, thrust::min(Fmax, F_cond));
		// 	const value_type G = thrust::max(0.0, thrust::min(Gmax, G_cond));
		// 	const value_type H = thrust::max(0.0, thrust::min(Hmax, H_cond));

		// 	const value_type lapl_u = u_top + u_bot + u_left + u_right - 4 * u;
		// 	const value_type lapl_v = v_top + v_bot + v_left + v_right - 4 * v;
		// 	const value_type lapl_w = w_top + w_bot + w_left + w_right - 4 * w;

		// 	// the dynamical equation
		// 	thrust::get<0>(thrust::get<8>(t)) = F - cu * u + Du * lapl_u * P_sin_theta;
		// 	thrust::get<1>(thrust::get<8>(t)) = G - cv * v + Dv * lapl_v * P_sin_theta;
		// 	thrust::get<2>(thrust::get<8>(t)) = H - cw * w + Dw * lapl_w * P_sin_theta;
		// }
	// };

	rd_dynamics(
		const state_type &init,
		const size_t &Nx_in, const size_t &Ny_in,
		const double &cu_in, const double &cv_in, const double &cw_in,
		const double &c1_in, const double &c2_in, const double &c3_in, const double &c4_in, const double &c5_in, const double &c6_in, const double &c7_in, const double &c8_in, const double &c9_in,
		const double &Du_in, const double &Dv_in, const double &Dw_in,
		const double &Fmax_in, const double &Gmax_in, const double &Hmax_in
	); //:
	// 	state_n (init) , N ( Nx_in * Ny_in ), Nx( Nx_in ), Ny ( Ny_in ),
	// 	cu(cu_in), cv(cv_in), cw(cw_in),
	// 	c1(c1_in), c2(c2_in), c3(c3_in), c4(c4_in), c5(c5_in), c6(c6_in), c7(c7_in), c8(c8_in), c9(c9_in),
	// 	Du(Du_in), Dv(Dv_in), Dw(Dw_in),
	// 	Fmax(Fmax_in), Gmax(Gmax_in), Hmax(Hmax_in),
	// 	top( init.size() ), bot( init.size() ), left( init.size() ), right( init.size() )
	// {
	// 	// define neighbours
	// 	thrust::counting_iterator<size_t> counter( 0 );

	// 	// top neighbours
	// 	thrust::copy( counter , counter+(N-Nx) , top.begin()+Nx );
	// 	thrust::copy( counter+(N-Nx), counter+N , top.begin() );

	// 	// bottom neighbours
	// 	thrust::copy( counter+Nx , counter+N , bot.begin() );
	// 	thrust::copy( counter, counter+Nx , bot.begin()+N-Nx );

	// 	// left neighbours
	// 	thrust::copy( counter , counter+N-1 , left.begin()+1 );

	// 	// right neighbours
	// 	thrust::copy( counter+1 , counter+N , right.begin() );

	// 	// adjust left and right neighbours on sides
	// 	for (int i = 0; i < Ny; ++i)
	// 	{
	// 		left[i * Nx] = i * Nx - 1;
	// 		right[(i + 1) * Nx - 1] = i * Nx;
	// 	}
	// }

	void operator() ( const state_type &x , state_type &dxdt , const value_type dt );

private:

	const state_type &state_n;
	const size_t N, Nx, Ny;
	const double cu, cv, cw;
	const double c1, c2, c3, c4, c5, c6, c7, c8, c9;
	const double Du, Dv, Dw;
	const double Fmax, Gmax, Hmax;
	index_vector_type top, bot, left, right;
};

void simulate_rd(Parameters &params);
