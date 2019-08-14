/*
 Copyright 2011-2013 Mario Mulansky
 Copyright 2011 Karsten Ahnert
 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

/*
 * This example shows how to use odeint on CUDA devices with thrust.
 * Note that we require at least Version 3.2 of the nVidia CUDA SDK
 * and the thrust library should be installed in the CUDA include
 * folder.
 *
 * As example we use a chain of phase oscillators with nearest neighbour
 * coupling, as described in:
 *
 * Avis H. Cohen, Philip J. Holmes and Richard H. Rand:
 * JOURNAL OF MATHEMATICAL BIOLOGY Volume 13, Number 3, 345-369,
 *
 */

#include <iostream>
#include <cmath>

#include <thrust/device_vector.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

using namespace std;

using namespace boost::numeric::odeint;


//change this to float if your device does not support double computation
typedef double value_type;


//[ thrust_phase_chain_system
//change this to host_vector< ... > if you want to run on CPU

// typedef thrust::tuple< value_type, value_type > tuple_type;
// typedef thrust::device_vector< tuple_type > state_type;

typedef thrust::device_vector< value_type > state_type;
typedef thrust::device_vector< size_t > index_vector_type;

//typedef thrust::host_vector< value_type > state_type;
//typedef thrust::host_vector< size_t > index_vector_type;

//<-
/*
 * This implements the rhs of the dynamical equation:
 * \phi'_0 = \omega_0 + sin( \phi_1 - \phi_0 )
 * \phi'_i  = \omega_i + sin( \phi_i+1 - \phi_i ) + sin( \phi_i - \phi_i-1 )
 * \phi'_N-1 = \omega_N-1 + sin( \phi_N-1 - \phi_N-2 )
 */
//->
class phase_oscillators
{

public:

    struct sys_functor
    {
        template< class Tuple >
        __host__ __device__
        void operator()( Tuple t )  // this functor works on tuples of values
        {
            // // first, unpack the tuple into value, neighbors and omega
            // const value_type phi = thrust::get<0>(t);
            // const value_type phi_left = thrust::get<2>(thrust::get<1>(t));  // left neighbor
            // const value_type phi_right = thrust::get<3>(thrust::get<1>(t)); // right neighbor
            // const value_type omega = thrust::get<7>(t);
            // // the dynamical equation
            // thrust::get<8>(t) = omega + sin( phi_right - phi ) + sin( phi - phi_left );

            // get component
            // const size_t component = thrust::get<9>(t);

            // get current values
            const value_type u = thrust::get<0>(thrust::get<0>(t));
            const value_type v = thrust::get<1>(thrust::get<0>(t));
            const value_type w = thrust::get<2>(thrust::get<0>(t));

            // get P sin(theta)
            const value_type P_sin_theta = thrust::get<3>(thrust::get<0>(t));

            // get neighbors
            const value_type u_top = thrust::get<0>(thrust::get<1>(t));  // top neighbor for u
            const value_type v_top = thrust::get<1>(thrust::get<1>(t));  // top neighbor for v
            const value_type w_top = thrust::get<2>(thrust::get<1>(t));  // top neighbor for w
            const value_type u_bot = thrust::get<0>(thrust::get<2>(t));  // bottom neighbor for u
            const value_type v_bot = thrust::get<1>(thrust::get<2>(t));  // bottom neighbor for v
            const value_type w_bot = thrust::get<2>(thrust::get<2>(t));  // bottom neighbor for w
            const value_type u_left = thrust::get<0>(thrust::get<3>(t));  // left neighbor for u
            const value_type v_left = thrust::get<1>(thrust::get<3>(t));  // left neighbor for v
            const value_type w_left = thrust::get<2>(thrust::get<3>(t));  // left neighbor for w
            const value_type u_right = thrust::get<0>(thrust::get<4>(t)); // right neighbor for u
            const value_type v_right = thrust::get<1>(thrust::get<4>(t)); // right neighbor for v
            const value_type w_right = thrust::get<2>(thrust::get<4>(t)); // right neighbor for w

            // get cu, cv, cw
            const value_type cu = thrust::get<0>(thrust::get<5>(t));
            const value_type cv = thrust::get<1>(thrust::get<5>(t));
            const value_type cw = thrust::get<2>(thrust::get<5>(t));

            // get Du, Dv, Dw
            const value_type Du = thrust::get<3>(thrust::get<5>(t));
            const value_type Dv = thrust::get<4>(thrust::get<5>(t));
            const value_type Dw = thrust::get<5>(thrust::get<5>(t));

            // get c1, c2, ..., c9
            const value_type c1 = thrust::get<0>(thrust::get<6>(t));
            const value_type c2 = thrust::get<1>(thrust::get<6>(t));
            const value_type c3 = thrust::get<2>(thrust::get<6>(t));
            const value_type c4 = thrust::get<3>(thrust::get<6>(t));
            const value_type c5 = thrust::get<4>(thrust::get<6>(t));
            const value_type c6 = thrust::get<5>(thrust::get<6>(t));
            const value_type c7 = thrust::get<6>(thrust::get<6>(t));
            const value_type c8 = thrust::get<7>(thrust::get<6>(t));
            const value_type c9 = thrust::get<8>(thrust::get<6>(t));

            // get Fmax, Gmax, Hmax
            const value_type Fmax = thrust::get<0>(thrust::get<7>(t));
            const value_type Gmax = thrust::get<1>(thrust::get<7>(t));
            const value_type Hmax = thrust::get<2>(thrust::get<7>(t));

            // get P
            // const value_type P = thrust::get<3>(thrust::get<7>(t));

            // compute each term for each component
            const value_type F_cond = c1 * v + c2 * w + c3;
            const value_type G_cond = c4 * u + c5 * w + c6;
            const value_type H_cond = c7 * u + c8 * v + c9;

            const value_type F = thrust::max(0.0, thrust::min(Fmax, F_cond));
            const value_type G = thrust::max(0.0, thrust::min(Gmax, G_cond));
            const value_type H = thrust::max(0.0, thrust::min(Hmax, H_cond));

            const value_type lapl_u = u_top + u_bot + u_left + u_right - 4 * u;
            const value_type lapl_v = v_top + v_bot + v_left + v_right - 4 * v;
            const value_type lapl_w = w_top + w_bot + w_left + w_right - 4 * w;

            // the dynamical equation
            thrust::get<0>(thrust::get<8>(t)) = F - cu * u + Du * lapl_u * P_sin_theta;
            thrust::get<1>(thrust::get<8>(t)) = G - cv * v + Dv * lapl_v * P_sin_theta;
            thrust::get<2>(thrust::get<8>(t)) = H - cw * w + Dw * lapl_w * P_sin_theta;
        }
    };

    phase_oscillators(
        const state_type &omega,
        const size_t &Nx_in, const size_t &Ny_in,
        const double &cu_in, const double &cv_in, const double &cw_in,
        const double &c1_in, const double &c2_in, const double &c3_in, const double &c4_in, const double &c5_in, const double &c6_in, const double &c7_in, const double &c8_in, const double &c9_in,
        const double &Du_in, const double &Dv_in, const double &Dw_in,
        const double &Fmax_in, const double &Gmax_in, const double &Hmax_in //,
        // const double &P_in
    ): m_omega( omega ) , m_N( omega.size() ) , m_prev( omega.size() ) , m_next( omega.size() ),
        state_n (omega) , N ( Nx_in * Ny_in ), Nx( Nx_in ), Ny ( Ny_in ),
        cu(cu_in), cv(cv_in), cw(cw_in),
        c1(c1_in), c2(c2_in), c3(c3_in), c4(c4_in), c5(c5_in), c6(c6_in), c7(c7_in), c8(c8_in), c9(c9_in),
        Du(Du_in), Dv(Dv_in), Dw(Dw_in),
        Fmax(Fmax_in), Gmax(Gmax_in), Hmax(Hmax_in),
        // P(P_in),
        top( omega.size() ), bot( omega.size() ), left( omega.size() ), right( omega.size() )//,
        // component( 3 * Nx_in * Ny_in )
    {
        // build indices pointing to left and right neighbours
        thrust::counting_iterator<size_t> c( 0 );
        thrust::copy( c , c+m_N-1 , m_prev.begin()+1 );
        m_prev[0] = 0; // m_prev = { 0 , 0 , 1 , 2 , 3 , ... , N-1 }

        thrust::copy( c+1 , c+m_N , m_next.begin() );
        m_next[m_N-1] = m_N-1; // m_next = { 1 , 2 , 3 , ... , N-1 , N-1 }

        // define neighbours
        thrust::counting_iterator<size_t> counter( 0 );

        // top neighbours
        thrust::copy( counter , counter+(N-Nx) , top.begin()+Nx );
        thrust::copy( counter+(N-Nx), counter+N , top.begin() );

        // bottom neighbours
        thrust::copy( counter+Nx , counter+N , bot.begin() );
        thrust::copy( counter, counter+Nx , bot.begin()+N-Nx );

        // left neighbours
        thrust::copy( counter , counter+N-1 , left.begin()+1 );

        // right neighbours
        thrust::copy( counter+1 , counter+N , right.begin() );

        // adjust left and right neighbours on sides
        for (int i = 0; i < Ny; ++i)
        {
            left[i * Nx] = i * Nx - 1;
            right[(i + 1) * Nx - 1] = i * Nx;
        }

        // set component values
        // thrust::fill_n(component.begin(), N, 0);
        // thrust::fill_n(component.begin()+N, N, 1);
        // thrust::fill_n(component.begin()+2*N, N, 2);
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

    const state_type &m_omega;
    const size_t m_N;
    index_vector_type m_prev;
    index_vector_type m_next;

    const state_type &state_n;
    const size_t N, Nx, Ny;
    const double cu, cv, cw;
    const double c1, c2, c3, c4, c5, c6, c7, c8, c9;
    const double Du, Dv, Dw;
    const double Fmax, Gmax, Hmax;
    // const double P;
    index_vector_type top, bot, left, right;
    // index_vector_type component;
};
//]

const size_t N = 8192;
const size_t Nx = 128, Ny = 64;
const value_type pi = 3.1415926535897932384626433832795029;
const value_type epsilon = 6.0 / ( N * N ); // should be < 8/N^2 to see phase locking
const value_type dt = 0.015;

int main( int arc , char* argv[] )
{
    // define constants
    const double cu=0.02, cv=0.025, cw=0.06;
    const double c1=-0.04, c2=-0.056, c3=0.382, c4=-0.05, c5=0, c6=0.25, c7=0.016, c8=-0.03, c9=0.24;
    const double Du=1.125, Dv=1.125, Dw=13.5;
    const double Fmax=0.5, Gmax=0.5, Hmax=0.5;
    const double P=0.00889;

    // create vcetors of data: all variables are concatenated into one vector for simplicity
    // create initial conditions and omegas on host:
    vector< value_type > x_host( 4 * N, 0 );
    vector< value_type > omega_host( 4 * N, 0 );
    for( size_t i=0 ; i<(3 * N) ; ++i )
    {
        x_host[i] = 2.0 * pi * drand48();
        omega_host[i] = ( 4 * N - i ) * epsilon; // decreasing frequencies
    }
    for( size_t i=3 * N ; i<(4 * N) ; ++i )
    {
        x_host[i] = 1 + P;
        omega_host[i] = 1 + P;
    }

    // copy to device
    state_type x = x_host;
    state_type omega = omega_host;

    // create stepper
    runge_kutta4< state_type , value_type , state_type , value_type > stepper;

    // create phase oscillator system function
    phase_oscillators sys(
        omega,
        Nx, Ny,
        cu, cv, cw,
        c1, c2, c3, c4, c5, c6, c7, c8, c9,
        Du, Dv, Dw,
        Fmax, Gmax, Hmax
    );

    // integrate
    integrate_const( stepper , sys , x , 0.0 , 10.0 , dt );

    // print some result
    thrust::copy( x.begin() , x.begin() + 10 , std::ostream_iterator< value_type >( std::cout , "\n" ) );
    std::cout << std::endl;
    //]
}

/*int main(int argc, const char *argv[])
{
    // Read parameters

    // Define variables
    int nx, ny;

    // Create initial state

    // Compute evolution

    // Export results
}*/
