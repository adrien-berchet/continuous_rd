/**
 * \file hex_2d_lattice.hpp
 *
 * \brief Provides a class to create and manipulate hexagonal 2D lattice objects.
 *
*/
#ifndef HEX_2D_LATTICE_HPP
#define HEX_2D_LATTICE_HPP
#include <fstream>
#include <memory>
#include <cmath>

#include <boost/math/constants/constants.hpp>

namespace generic
{
	typedef thrust::host_vector< int > host_int_vector_type;

	const double PI = boost::math::constants::pi<double>();

	namespace _internals
	{
		struct Point {
			double x, y;

			Point(const double &x_in, const double &y_in): x(x_in), y(y_in) {}
		};

		struct Segment {
			Point p1, p2;

			Segment(const Point &p1_in, const Point &p2_in): p1(p1_in), p2(p2_in) {}
		};
	}

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
			// const double &side_in, const double &x_length_in, const double &y_length_in
		):
			side(side_in),
			Nq(Nq_in), Nr(Nr_in), N(Nq * Nr),
			hex_width(2. * side), hex_hor_space(hex_width * 3. / 4.), hex_vert_space(side * sqrt(3.)),
			x_length(Nq * hex_hor_space), y_length(Nr * hex_vert_space),
			x(generate_x(Nq, Nr)), y(generate_y(Nq, Nr))
		{}

		inline
		size_t eucl_to_hex_ind(const double &x, const double &y) const
		{
			return 0;
		}

		inline
		std::tuple<double, double> hex_coords_to_eucl(const double &q, const double &r) const
		{
			const double shift = (std::fmod(std::fmod(q, Nq), 2) == 0 ? 0. : 0.5 * hex_vert_space);
			return std::make_tuple(
				(q / Nq) * hex_hor_space + shift,
				r * hex_vert_space
			);
		}

		inline
		std::tuple<size_t, size_t> eucl_to_hex_coords(const double &x, const double &y) const
		{
			    double q = ( (2. / 3.) * x ) / side;
			    double r = (1. / 3.) * ( -x  +  (std::sqrt(3.) * y) ) / side;

				// Round coordinates
				int xi = int(std::round(q));
				int zi = int(std::round(r));
				int yi = int(std::round(-q - r));
				const double x_diff = std::abs(double(xi) - q);
				const double z_diff = std::abs(double(zi) - r);
				const double y_diff = std::abs(double(yi) - (-q - r));
				if (x_diff > z_diff && x_diff > y_diff)
				{
					xi = -zi - yi;
				}
				else if (y_diff > z_diff)
				{
					yi = -xi - zi;
				}
				else
				{
					zi = -xi - yi;
				}

				const size_t col = xi >= 0 ? xi : Nq - 1;
				int row = zi + int((xi - (xi & 1)) / 2.);

				if(row < 0)
				{
					row = Nr - 1;
				}

				return std::make_tuple(col, size_t(row));
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
		 * \brief Compute the angle between two segments. Return -10 if they do not intersect.
		 * \param s1 First segment.
		 * \param s2 Second segment.
		*/
		double intersection_sin_angle(
			const std::vector<double> &s1,
			const std::vector<double> &s2
		) const
		{
			_internals::Point s1_p1(s1[0], s1[1]);
			_internals::Point s1_p2(s1[2], s1[3]);
			_internals::Point s2_p1(s2[0], s2[1]);
			_internals::Point s2_p2(s2[2], s2[3]);
			_internals::Segment l1(s1_p1, s1_p2);
			_internals::Segment l2(s2_p1, s2_p2);

			// Four direction for two lines and points of other line
			int dir1 = direction(l1.p1, l1.p2, l2.p1);
			int dir2 = direction(l1.p1, l1.p2, l2.p2);
			int dir3 = direction(l2.p1, l2.p2, l1.p1);
			int dir4 = direction(l2.p1, l2.p2, l1.p2);

			int intersects = 0;
			if(dir1 != dir2 && dir3 != dir4)
				intersects = 1; // They are intersecting

			if(dir1==0 && onLine(l1, l2.p1)) // When p2 of line2 is on the line1
				intersects = 2;

			if(dir2==0 && onLine(l1, l2.p2)) // When p1 of line2 is on the line1
				intersects = 2;

			if(dir3==0 && onLine(l2, l1.p1)) // When p2 of line1 is on the line2
				intersects = 2;

			if(dir4==0 && onLine(l2, l1.p2)) // When p1 of line1 is on the line2
				intersects = 2;

			// Compute angle
			if(intersects == 2)
			{
				return 0;
			}
			else if(intersects == 1)
			{
				const double x1 = (s1_p2.x - s1_p1.x);
				const double x2 = (s2_p2.x - s2_p1.x);
				const double y1 = (s1_p2.y - s1_p1.y);
				const double y2 = (s2_p2.y - s2_p1.y);
				const double det = (x1 * y2 - y1 * x2);
				const double norm = std::sqrt(std::pow(x1, 2) + std::pow(y1, 2)) * std::sqrt(std::pow(x2, 2) + std::pow(y2, 2));
				return (norm != 0 ? det / norm : -10);
			}
			else
			{
				return -10;
			}
		}

	public: // Attributes
		const double side;
		const size_t Nq, Nr, N;
		const double hex_width;
		const double hex_hor_space;
		const double hex_vert_space;
		const double x_length, y_length;
		std::vector<double> x, y;

	private: // Methods

		std::vector<double> generate_x(const size_t &Nq, const size_t &Nr)
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

		std::vector<double> generate_y(const size_t &Nq, const size_t &Nr)
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


		bool onLine(_internals::Segment l1, _internals::Point p) const
		{
			//check whether p is on the line or not
		   if(p.x <= std::max(l1.p1.x, l1.p2.x) && p.x <= std::min(l1.p1.x, l1.p2.x) &&
		      (p.y <= std::max(l1.p1.y, l1.p2.y) && p.y <= std::min(l1.p1.y, l1.p2.y)))
		      return true;

		   return false;
		}

		int direction(_internals::Point a, _internals::Point b, _internals::Point c) const
		{
			int val = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);
			if (val == 0)
		    	return 0;     //colinear
		  	else if(val < 0)
		    	return 2;    //anti-clockwise direction
		    else
	    		return 1;    //clockwise direction
		}

	// private: // Attributes
	};

/**
 * \brief Export hexagonal lattice. Only used for validation.
*/
template<typename T>
void export_hex_lattice(const Hex2dLattice &hex, const T &hex_colors, const std::string &folder)
{
	BOOST_LOG_TRIVIAL(info) << "Export hexagonal lattice";

	// Create file
	generic::DatWriter data_file(folder + "/hex_lattice.dat");

	// Write header
	data_file.write_header("Hexagonal lattice", hex.Nq, hex.Nr, "x", "y", "color");

	// Write data
	for(
		auto i=thrust::make_zip_iterator(
			thrust::make_tuple(
				hex.x.begin(),
				hex.y.begin(),
				hex_colors.begin()
		) );
		i != thrust::make_zip_iterator(
			thrust::make_tuple(
				hex.x.end(),
				hex.y.end(),
				hex_colors.end()
		) );
		++i
	)
	{
		data_file.write_row(
			thrust::get<0>(*i),
			thrust::get<1>(*i),
			thrust::get<2>(*i)
		);
	}
}

} // namespace generic

#endif // HEX_2D_LATTICE_HPP
