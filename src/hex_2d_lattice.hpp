/**
 * \file hex_2d_lattice.hpp
 *
 * \brief Provides a class to create and manipulate hexagonal 2D lattice objects.
 *
*/
#ifndef HEX_2D_LATTICE_HPP
#define HEX_2D_LATTICE_HPP

#include <boost/math/constants/constants.hpp>

#include "dat_writer.hpp"
#include "logger.hpp"

namespace generic
{
	/** \brief Pi value. */
	const double PI = boost::math::constants::pi<double>();

	/**
	 * \ingroup generic
	 * \brief This is a class that provides methods to create and manipulate hexagonal 2D lattice objects with vertical layout.
	 * The implementation was inspired by https://www.redblobgames.com/grids/hexagons.
	*/
	class Hex2dLattice
	{
	public: // Methods

		/**
		 * \brief Constructor of the Hex2dLattice.
		 * \param side_in Length of a hexagonal side.
		 * \param Nq_in The size along the X axis of the lattice.
		 * \param Nr_in The size along the Y axis of the lattice.
		*/
		Hex2dLattice(
			const double &side_in, const size_t &Nq_in, const size_t &Nr_in
		):
			side(side_in),
			Nq(Nq_in), Nr(Nr_in), N(Nq * Nr),
			hex_width(2. * side), hex_hor_space(hex_width * 3. / 4.), hex_vert_space(side * sqrt(3.)),
			x_length(Nq * hex_hor_space), y_length(Nr * hex_vert_space),
			x(generate_x()), y(generate_y()), half_side(0.5 * side_in)
		{}

		/**
		 * \brief Convert hexagonal coordinates to euclidean coordinates.
		*/
		inline
		std::tuple<double, double> hex_coords_to_eucl(const double &q, const double &r) const
		{
			const double remainder = std::fmod(q, 2);
			const double shift = (remainder <= 1 ? remainder * 0.5 * hex_vert_space : (2. - remainder) * 0.5 * hex_vert_space);
			return std::make_tuple(
				q * hex_hor_space,
				r * hex_vert_space + shift
			);
		}

		/**
		 * \brief Convert euclidean coordinates to rounded hexagonal coordinates.
		 * \param x The X coordinate.
		 * \param y The Y coordinate.
		*/
		inline
		std::tuple<size_t, size_t> eucl_to_rounded_hex_coords(const double &x, const double &y) const
		{
			// Compute axial coordinates
				const double q = (2. * x ) / (3. * side);
				const double r = (1. / 3.) * ( -x  +  (std::sqrt(3.) * y) ) / side;

			// Round cube coordinates to find the nearest hexagon center
			int xi = int(std::round(q));
			int zi = int(std::round(r));
			int yi = int(std::round(-q - r));
			const double x_diff = std::abs(double(xi) - q);
			const double z_diff = std::abs(double(zi) - r);
			const double y_diff = std::abs(double(yi) - (-q - r));
			BOOST_LOG_TRIVIAL(trace) << "("<<x<<", "<<y<<") => (q, r) = ("<<q<<", "<<r<<") ; (xi, zi, yi) = ("<<xi<<", "<<zi<<", "<<yi<<") ; diff = ("<<x_diff<<", "<<z_diff<<", "<<y_diff<<")";

			// Adjust values if needed
			if (z_diff > x_diff && z_diff > y_diff)
			{
				zi = -xi - yi;
			}
			else if (x_diff > y_diff)
			{
				xi = -zi - yi;
			}
			else
			{
				yi = -xi - zi;
			}

			// Convert back to col and row coordinates
			const int col = xi;
			const int row = zi + size_t(xi / 2.);

			BOOST_LOG_TRIVIAL(trace) << " round(xi, zi, yi) = ("<<xi<<", "<<zi<<", "<<yi<<") ; border(xi, zi, yi) = ("<<xi<<", "<<zi<<", "<<yi<<")"<<" ; col, row = "<<col<<", "<<row<<std::endl;

			// Force positive values
			return std::make_tuple(
				col >= int(Nq) ? col - Nq : (col < 0 ? col + Nq : col),
				row >= int(Nr) ? row - Nr : (row < 0 ? row + Nr : row)
			);
		}

		/**
		 * \brief Get hexagonal coordinates from index.
		 * \param ind Index of the hexagon in the lattice.
		*/
		inline
		std::tuple<size_t, size_t> ind_to_hex_coords(const size_t &ind) const
		{
			return std::make_tuple(
				size_t(ind % Nq),
				size_t(ind / Nq)
			);
		}

		/**
		 * \brief Get index from hexagonal coordinates.
		 * \param coords Tuplus of coordinates.
		*/
		inline
		size_t hex_coords_to_ind(const std::tuple<size_t, size_t> &coords) const
		{
			return std::get<0>(coords) + Nq * std::get<1>(coords);
		}

		/**
		 * \brief Build side segments of a hexagon.
		 * \param ind Index of the hexagon in the lattice.
		*/
		std::vector<std::vector<double>> hex_segments(const size_t &ind) const
		{
			std::vector<std::vector<double>> segments;
			for (size_t i = 0; i < 6; ++i)
			{
				const double angle_rad = PI / 180. * (60. * i);
				const double angle_rad_next = PI / 180. * (60. * (i + 1));
				std::vector<double> side_segment{
					x[ind] + side * std::cos(angle_rad),
					y[ind] + side * std::sin(angle_rad),
					x[ind] + side * std::cos(angle_rad_next),
					y[ind] + side * std::sin(angle_rad_next)};
				segments.push_back(std::move(side_segment));
			}
			return segments;
		}

		/**
		* \brief Compute the P(xx') terms of a given point.
		* \param ind Index of the hexagon in which the point is located.
		* \param x X coordinate of the point.
		* \param y Y coordinate of the point.
		* \param epsilon Size of the buffer in which P(xx') < 1.
		* \param P The factor applied to the Laplacian in the buffer.
		*/
		std::tuple<double, double, double, double> Pxx_terms(const size_t &ind, const double &x, const double &y, const double &epsilon, const double &P) const
		{
			std::tuple<double, double, double, double> terms{1, 1, 1, 1};
			const double &hex_x = this->x[ind];
			const double &hex_y = this->y[ind];

			// Top and bot components
			if (x <= hex_x - half_side)
			{
				// Left part
				if (y >= sqrt3 * x - sqrt3 * (hex_x - side) + hex_y - epsilon)
				{
					std::get<0>(terms) = P * sin30;
				}
				if (y <= -sqrt3 * x + sqrt3 * (hex_x - side) + hex_y + epsilon)
				{
					std::get<1>(terms) = P * sin30;
				}
			}
			else if (x >= hex_x + half_side)
			{
				// Right part
				if (y >= -sqrt3 * x + sqrt3 * (hex_x + side) + hex_y - epsilon)
				{
					std::get<0>(terms) = P * sin30;
				}
				if (y <= sqrt3 * x - sqrt3 * (hex_x + side) + hex_y + epsilon)
				{
					std::get<1>(terms) = P * sin30;
				}
			}
			else
			{
				// Middle part
				if (y > hex_y + 0.5 * hex_vert_space - epsilon)
				{
					std::get<0>(terms) = P;
				}
				else if (y <= hex_y - 0.5 * hex_vert_space + epsilon)
				{
					std::get<1>(terms) = P;
				}
			}

			// Left and right components
			if (x <= hex_x - half_side + epsilon)
			{
				// Left part
				if ((y >= sqrt3 * x - sqrt3 * (hex_x - side) + hex_y - sqrt3 * epsilon) ||
					(y <= -sqrt3 * x + sqrt3 * (hex_x - side) + hex_y + sqrt3 * epsilon))
				{
					std::get<2>(terms) = P * sin60;
				}
			}
			else if (x >= hex_x + half_side - epsilon)
			{
				// Right part
				if ((y >= -sqrt3 * x + sqrt3 * (hex_x + side) + hex_y - sqrt3 * epsilon) ||
					(y <= sqrt3 * x - sqrt3 * (hex_x + side) + hex_y + sqrt3 * epsilon))
				{
					std::get<3>(terms) = P * sin60;
				}
			}

			return terms;
		}

	public: // Attributes
		/**@{ */
		const double side;
		const size_t Nq, Nr, N;
		const double hex_width;
		const double hex_hor_space;
		const double hex_vert_space;
		const double x_length, y_length;
		std::vector<double> x, y;
		/**@}*/

	private: // Methods

		/**
		 * \brief Build the X coordinates of the hexagons.
		*/
		std::vector<double> generate_x()
		{
			std::vector<double> x_tmp(Nq * Nr);

			std::vector<double> x_coords(Nq);
			for (size_t i = 0; i < Nq; ++i)
			{
				x_coords[i] = i * this->hex_hor_space;
			}

			// Create X coordinates
			for (size_t i = 0; i < Nr; ++i)
			{
				std::copy(x_coords.begin(), x_coords.end(), x_tmp.begin() + i * Nq);
			}

			return x_tmp;
		}

		/**
		 * \brief Build the Y coordinates of the hexagons.
		*/
		std::vector<double> generate_y()
		{
			std::vector<double> y_tmp(Nq * Nr);

			// Create Y coordinates
			for (size_t i = 0; i < Nq * Nr; ++i)
			{
				const double shift = ((i % Nq) % 2 == 0 ? 0 : 0.5 * hex_vert_space);
				const double y_i = size_t(i / Nq) * hex_vert_space + shift;
				y_tmp[i] = std::move(y_i);
			}

			return y_tmp;
		}

	private:
		/**@{ */
		const double half_side;
		const double sqrt3{std::sqrt(3.)};
		const double sin30{0.5};
		const double sin60{0.5 * std::sqrt(3.)};
		/**@}*/
	};

/**
 * \brief Export hexagonal lattice. Only used for validation.
*/
template<typename T>
void export_hex_lattice(const Hex2dLattice &hex, const T &hex_colors, const std::string &folder)
{
	std::string filename = folder + "/hex_lattice.dat";
	BOOST_LOG_TRIVIAL(info) << "Export hexagonal lattice to "<<filename;

	// Create file
	generic::DatWriter data_file(filename);

	// Write header
	data_file.write_header("Hexagonal lattice", hex.Nq, hex.Nr, "x", "y", "color");

	// Write data
	// Use raw iterators because of issues with Boost::test
	auto x = hex.x.begin();
	auto y = hex.y.begin();
	auto colors = hex_colors.begin();
	auto x_end = hex.x.end();
	auto y_end = hex.y.end();
	auto colors_end = hex_colors.end();
	for (; x != x_end, y != y_end, colors != colors_end; ++ x, ++y, ++colors)
	{
		data_file.write_row(
			*x,
			*y,
			*colors
		);
	}
}

} // namespace generic

#endif // HEX_2D_LATTICE_HPP
