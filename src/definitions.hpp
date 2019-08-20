/**
 * \file definitions.hpp
 *
 * \brief Provides some global typedefs, global constants and functions.
 *
*/
#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP


#include <thrust/device_vector.h>


/** Default values for black scales. */
const std::tuple<double, double, double> black(6.85, 0.09, 4.75);
/** Default values for green scales. */
const std::tuple<double, double, double> green(0.05, 5.35, 0.09);


/** Default value type. Change this to float if your device does not support double computation. */
typedef double value_type;


/** Default host vector type that store velue_type values. This type is used to store the states of dynamic systems on the host. */
typedef thrust::host_vector< value_type > host_state_type;

/** Default host vector type that store integer values. */
typedef thrust::host_vector< int > host_int_vector_type;

/** Default host vector type that store size_t values. */
typedef thrust::host_vector< size_t > host_index_vector_type;

/** Default device vector type that store value_type values. This type is used to store the states of dynamic systems on the device. */
typedef thrust::device_vector< value_type > state_type;

/** Default device vector type that store size_t values. */
typedef thrust::device_vector< size_t > index_vector_type;


/**
 * \brief Compute X coordinate
*/
inline
double get_x(const size_t &it, const size_t &Nx, const double &epsilon)
{
	return (it % Nx) * epsilon;
}


/**
 * \brief Compute X coordinate
*/
inline
double get_y(const size_t &it, const size_t &Nx, const double &epsilon)
{
	return int(it / Nx) * epsilon;
}

#endif // DEFINITIONS_HPP
