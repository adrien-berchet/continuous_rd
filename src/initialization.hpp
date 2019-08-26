/**
 * \file initialization.hpp
 *
 * \brief Provides tools to initialize the dynamic system.
 *
*/
#ifndef INITIALIZATION_HPP
#define INITIALIZATION_HPP

#include <random>

#include "class_parameters.hpp"
#include "hex_2d_lattice.hpp"

/**
 * \brief Random initialization
*/
template<typename T>
void random_hex_init(T &state, const Parameters &params, T &Pxx_top, T &Pxx_bot, T &Pxx_left, T &Pxx_right)
{
	const size_t &Nx = params.Nx, &Ny = params.Ny;
	const size_t N = Nx * Ny;

	const double &epsilon = params.epsilon;
	const double &S = params.S;
	const double &P = params.P;

	BOOST_LOG_TRIVIAL(info) << "Random initialization:";
	BOOST_LOG_TRIVIAL(info) << "\t Nx = "<<Nx;
	BOOST_LOG_TRIVIAL(info) << "\t Ny = "<<Ny;
	BOOST_LOG_TRIVIAL(info) << "\t epsilon = "<<epsilon;
	BOOST_LOG_TRIVIAL(info) << "\t S = "<<S;
	BOOST_LOG_TRIVIAL(info) << "\t P = "<<P;


	// Generate hexagonal lattice
	const double hex_width(2. * S), hex_hor_space(hex_width * 3. / 4.), hex_vert_space(S * sqrt(3.));
	const size_t Nq(size_t(1 + (Nx * epsilon) / hex_hor_space));
	const size_t Nr(size_t(1 + (Ny * epsilon) / hex_vert_space));

	generic::Hex2dLattice hex(S, Nq, Nr);

	BOOST_LOG_TRIVIAL(debug) << "\t Hexagonal lattice generation";
	BOOST_LOG_TRIVIAL(debug) << "\t\t hex.side="<<hex.side;
	BOOST_LOG_TRIVIAL(debug) << "\t\t hex.x_length="<<hex.x_length;
	BOOST_LOG_TRIVIAL(debug) << "\t\t hex.y_length="<<hex.y_length;
	BOOST_LOG_TRIVIAL(debug) << "\t\t hex.hex_width="<<hex.hex_width;
	BOOST_LOG_TRIVIAL(debug) << "\t\t hex.hex_hor_space="<<hex.hex_hor_space;
	BOOST_LOG_TRIVIAL(debug) << "\t\t hex.hex_vert_space="<<hex.hex_vert_space;
	BOOST_LOG_TRIVIAL(debug) << "\t\t hex.Nq="<<hex.Nq;
	BOOST_LOG_TRIVIAL(debug) << "\t\t hex.Nr="<<hex.Nr;

	std::vector<int> hex_colors(hex.N, 0);

    // Generate random colors
	std::random_device rd;
    std::mt19937 random_engine(rd());
    std::uniform_int_distribution<> distribution(0, 1);
	std::generate(hex_colors.begin(), hex_colors.end(), [&distribution, &random_engine]() { return distribution(random_engine); });

	BOOST_LOG_TRIVIAL(debug) << "\t Generate component values and Laplacian correction terms";

	// Initialize P * sin(theta) terms
	// Set P sin(theta) = epsilon^{-2} everywhere for all direction
	std::fill(Pxx_top.begin(), Pxx_top.end(), std::pow(epsilon, -2.0));
	std::fill(Pxx_bot.begin(), Pxx_bot.end(), std::pow(epsilon, -2.0));
	std::fill(Pxx_left.begin(), Pxx_left.end(), std::pow(epsilon, -2.0));
	std::fill(Pxx_right.begin(), Pxx_right.end(), std::pow(epsilon, -2.0));

    std::uniform_real_distribution<> perturbation(-0.1, 0.1);
	thrust::counting_iterator<size_t> num(0);
	for(
		auto i=thrust::make_zip_iterator(thrust::make_tuple(
			state.begin(),
			state.begin() + N,
			state.begin() + 2 * N,
			Pxx_top.begin(),
			Pxx_bot.begin(),
			Pxx_left.begin(),
			Pxx_right.begin(),
			num
		) );
		i != thrust::make_zip_iterator(thrust::make_tuple(
			state.begin() + N,
			state.begin() + 2 * N,
			state.begin() + 3 * N,
			Pxx_top.end(),
			Pxx_bot.end(),
			Pxx_left.end(),
			Pxx_right.end(),
			num + N
		) );
		++i
	)
	{
		// Compute coordinates of the current element
		const double x = get_x(thrust::get<7>(*i), Nx, epsilon);
		const double y = get_y(thrust::get<7>(*i), Nx, epsilon);

		// Find to which hexagon it belongs
		size_t hex_ind = hex.hex_coords_to_ind(hex.eucl_to_rounded_hex_coords(x, y));

		// Find the coordinates and color of this hexagon
		const int &color_hex = hex_colors[hex_ind];

		// Compute values of u, v, w
		const double u_init = (color_hex == 0 ? std::get<0>(black) : std::get<0>(green)) + perturbation(random_engine);
		const double v_init = (color_hex == 0 ? std::get<1>(black) : std::get<1>(green)) + perturbation(random_engine);
		const double w_init = (color_hex == 0 ? std::get<2>(black) : std::get<2>(green)) + perturbation(random_engine);

		// Set values to u, v, w
		thrust::get<0>(*i) = u_init;
		thrust::get<1>(*i) = v_init;
		thrust::get<2>(*i) = w_init;

		// Compute the Laplacian correction terms
		std::tuple<double, double, double, double> pxx_terms = hex.Pxx_terms(hex_ind, x, y, epsilon, P);
		thrust::get<3>(*i) *= std::get<0>(pxx_terms);
		thrust::get<4>(*i) *= std::get<1>(pxx_terms);
		thrust::get<5>(*i) *= std::get<2>(pxx_terms);
		thrust::get<6>(*i) *= std::get<3>(pxx_terms);
	}

	// Export hexagonal lattice if requested
	if(params.export_hex_lattice)
	{
		generic::export_hex_lattice(hex, hex_colors, params.result_folder);
	}
}

/**
 * \brief Gaussian initialization. Only used for validation.
*/
template<typename T>
void gauss_init(T &state, const Parameters &params, T &Pxx_top, T &Pxx_bot, T &Pxx_left, T &Pxx_right)
{
	const size_t &Nx = params.Nx, &Ny = params.Ny;
	const size_t N = Nx * Ny;

	const double &epsilon = params.epsilon;
	const double &sigma = params.gauss_std;

	const double x_center = (Nx / 2.0) * epsilon;
	const double y_center = (Ny / 2.0) * epsilon;

	BOOST_LOG_TRIVIAL(info) << "Gaussian initialization:";
	BOOST_LOG_TRIVIAL(info) << "\t Nx = "<<Nx;
	BOOST_LOG_TRIVIAL(info) << "\t Ny = "<<Ny;
	BOOST_LOG_TRIVIAL(info) << "\t x_center = "<<x_center;
	BOOST_LOG_TRIVIAL(info) << "\t y_center = "<<y_center;
	BOOST_LOG_TRIVIAL(info) << "\t epsilon = "<<epsilon;
	BOOST_LOG_TRIVIAL(info) << "\t std = "<<sigma;

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
		const double x = get_x(num, Nx, epsilon);
		const double y = get_y(num, Nx, epsilon);
		const double r = std::sqrt(std::pow(x - x_center, 2.) + std::pow(y - y_center, 2.));
		const double C = std::exp(-std::pow(r, 2.) / (2. * std::pow(sigma, 2.)));
		thrust::get<0>(*i) = C;
		thrust::get<1>(*i) = C;
		thrust::get<2>(*i) = C;
		++num;
	}

	// Set P sin(theta) = 1 everywhere for all direction
	std::fill(Pxx_top.begin(), Pxx_top.end(), 1.0);
	std::fill(Pxx_bot.begin(), Pxx_bot.end(), 1.0);
	std::fill(Pxx_left.begin(), Pxx_left.end(), 1.0);
	std::fill(Pxx_right.begin(), Pxx_right.end(), 1.0);
}

#endif // INITIALIZATION_HPP
