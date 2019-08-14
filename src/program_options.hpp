/**
 * \file program_options.hpp
 *
 * \brief Provides functions to load command line parameters
 * or to read a file parameter.
 *
*/
#ifndef PROGRAM_OPTIONS_HPP
#define PROGRAM_OPTIONS_HPP
#include <iostream>
#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace generic
{

	/**
	 * \ingroup Generic
	 * \brief Transforms any standard type to string
	*/
	template<class T>
	inline
	std::string new_string(T &a)
	{
		return static_cast<std::ostringstream*>( &(std::ostringstream() << a) )->str();
	}

	/**
	 * \ingroup Generic
	 * \brief This is a container to store and process a vector of keys.
	*/
	class MultiKey
	{
	public:
		/**
		 * \brief Vector of keys.
		*/
		std::vector<std::string> keys;

		/**
		 * \brief Constructor of the MultiKey.
		 * \param a String added to the vector of keys.
		*/
		MultiKey(const std::string &a)
		{
			keys.push_back(a);
		}

		/**
		 * \brief Constructor of the MultiKey.
		 * \param a Vector of keys.
		*/
		MultiKey(const std::vector<std::string> &a): keys(a){}

		/**
		 * \brief Constructor of the MultiKey.
		 * \param a Vector of keys.
		*/
		MultiKey(std::vector<std::string> &&a)
		{
			keys.swap(a);
		}

		/**
		 * \brief Gets the number of keys.
		 * \return Size of the vector of keys.
		*/
		size_t size() const
		{
			return this->keys.size();
		}

		/**
		 * \brief Adds a key to the vector of keys.
		 * \param a New key to add.
		*/
		void push_back(const std::string &a)
		{
			keys.push_back(a);
		}

		/**
		 * \brief Gets the first key.
		 * \return A string containing the first key.
		*/
		std::string front() const
		{
			return keys.front();
		}

		/**
		 * \brief Gets the last key.
		 * \return A string containing the last key.
		*/
		std::string back() const
		{
			return keys.back();
		}

		/**
		 * \brief Compares a key to all keys contained in the MultiKey.
		 * \param str The string to compare.
		 * \return An int that is 0 when the key is equal to a key stored in the MultiKey.
		*/
		int compare(const std::string &str) const
		{
			int ok=0;
			std::string tmp1=str;
			std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::tolower);
			for (size_t i = 0; i < keys.size(); ++i)
			{
				std::string tmp2=keys[i];
				std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::tolower);
				ok=tmp2.compare(tmp1);
				if (ok==0)
				{
					return 0;
				}
			}
			return ok;
		}

		/**
		 * \brief Compares a key to all keys contained in the MultiKey.
		 * \param pos The first character to compare.
		 * \param len The number of characters to compare.
		 * \param str The string to compare.
		 * \return An int that is 0 when the key is equal to a key stored in the MultiKey.
		*/
		int compare(const size_t &pos, const size_t &len , const std::string &str) const
		{
			int ok=0;
			std::string tmp1=str;
			std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::tolower);
			for (size_t i = 0; i < keys.size(); ++i)
			{
				std::string tmp2=keys[i].substr(pos, len);
				std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::tolower);
				ok=tmp2.compare(tmp1);
				if (ok==0)
				{
					return 0;
				}
			}
			return ok;
		}

		/**
		 * \brief Compares a key to all keys contained in the MultiKey.
		 * \param str The string to compare.
		 * \return An int that is 0 when the key is equal to a key stored in the MultiKey.
		*/
		int compare(const char* str) const
		{
			int ok=0;
			std::string tmp1=str;
			std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::tolower);
			for (size_t i = 0; i < keys.size(); ++i)
			{
				std::string tmp2=keys[i];
				std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::tolower);
				ok=tmp2.compare(tmp1);
				if (ok==0)
				{
					return 0;
				}
			}
			return ok;
		}

		/**
		 * \brief Compares a key to all keys contained in the MultiKey.
		 * \param pos The first character to compare.
		 * \param len The number of characters to compare.
		 * \param str The string to compare.
		 * \return An int that is 0 when the key is equal to a key stored in the MultiKey.
		*/
		int compare(const size_t &pos, const size_t &len , const char* str) const
		{
			int ok=0;
			std::string tmp1=str;
			std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::tolower);
			for (size_t i = 0; i < keys.size(); ++i)
			{
				std::string tmp2=keys[i].substr(pos, len);
				std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::tolower);
				ok=tmp2.compare(tmp1);
				if (ok==0)
				{
					return 0;
				}
			}
			return ok;
		}

		/**
		 * \brief Adds a key to the MultiKey.
		 * \param a The key to add.
		*/
		void operator[](const std::string &a)
		{
			keys.push_back(a);
		}

		/**
		 * \brief Gets a key from the MultiKey.
		 * \param a The index of the key to get.
		 * \return A string containing the key.
		*/
		std::string operator[](const size_t &a) const
		{
			return keys[a];
		}

		/**
		 * \brief Checks if a key is contained in the MultiKey.
		 * \param a The key to look for.
		 * \return True or False.
		*/
		bool operator==(const std::string &a) const
		{
			std::string tmp1=a.substr(0,a.find_first_of(" "));
			std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::tolower);
			for (size_t i = 0; i < keys.size(); ++i)
			{
				std::string tmp2=keys[i].substr(0,keys[i].find_first_of(" "));
				std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::tolower);
				if (tmp1==tmp2)
				{
					return true;
				}
			}
			return false;
		}

		/**
		 * \brief Checks if a MultiKey is equal to the current MultiKey.
		 * \param a The MultiKey to check.
		 * \return True or False.
		*/
		bool operator==(const MultiKey &a) const
		{
			for (size_t i = 0; i < keys.size(); ++i)
			{
				std::string tmp1=keys[i].substr(0,keys[i].find_first_of(" "));
				std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::tolower);
				for (size_t j = 0; j < a.keys.size(); ++j)
				{
					std::string tmp2=a.keys[j].substr(0,a.keys[j].find_first_of(" "));
					std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::tolower);
					if (tmp1==tmp2)
					{
						return true;
					}
				}
			}
			return false;
		}

		/**
		 * \brief Structure derived from a std::binary_function used to compare two strings case insensitively.
		*/
		struct nocase_compare : public std::binary_function<unsigned char,unsigned char,bool>
		{
			/**
			 * \brief Operator used to compare two strings case insensitively.
			 * \param c1 The first string to compare.
			 * \param c2 The second string to compare.
			 * \return True or False.
			*/
			bool operator() (std::string c1, std::string c2) const
			{
				std::transform(c1.begin(), c1.end(), c1.begin(), ::tolower);
				std::transform(c2.begin(), c2.end(), c2.begin(), ::tolower);
				return c1 < c2;
			}
		};

		/**
		 * \brief Operator used to sort the MultiKey(s).
		 * \param a The string to compare.
		 * \return True or False.
		*/
		bool operator<(const std::string &a) const
		{
			std::string tmp2=a.substr(0,a.find_first_of(" "));
			std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::tolower);

			for (size_t i = 0; i < keys.size(); ++i)
			{
				std::string tmp1=keys[i].substr(0,keys[i].find_first_of(" "));
				std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::tolower);
				if (tmp1==tmp2)
					return false;
			}

			std::vector<std::string> tmp=keys;
			std::sort(tmp.begin(), tmp.end());

			std::string tmp1=tmp.front().substr(0,tmp.front().find_first_of(" "));
			std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::tolower);
			return tmp1 < tmp2;
		}

		/**
		 * \brief Operator used to sort the MultiKey(s).
		 * \param a The MultiKey to compare.
		 * \return True or False.
		*/
		bool operator<(const MultiKey &a) const
		{
			for (size_t i = 0; i < keys.size(); ++i)
			{
				std::string tmp1=keys[i].substr(0,keys[i].find_first_of(" "));
				std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::tolower);
				for (size_t j = 0; j < a.keys.size(); ++j)
				{
					std::string tmp2=a.keys[j].substr(0,a.keys[j].find_first_of(" "));
					std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::tolower);
					if (tmp1==tmp2)
						return false;
				}
			}

			std::string tmp1=keys.front().substr(0,keys.front().find_first_of(" "));
			std::string tmp2=a.keys.front().substr(0,a.keys.front().find_first_of(" "));
			std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::tolower);
			std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::tolower);
			return tmp1 < tmp2;
		}
	};

	/**
	 * \ingroup Generic
	 * \brief This is a container to store and process parameters.
	 * \deprecated This is quite deprecated now because Boost::program_options can replace it.
	*/
	class ProgramOptions
	{
	public:


		/**
		 * \brief Constructor of the ProgramOptions object.
		*/
		ProgramOptions()
		{
			code_name_="ProgramOptions";
			description_="Test ProgramOptions class.";
			synopsis_=code_name_+" -i [FILE] -o [FILE] [OPTIONS]...";
			initialize();
		};

		/**
		 * \brief Constructor of the ProgramOptions object.
		 * \param name Code name.
		 * \param desc Code description.
		 * \param syn Code synopsis.
		*/
		ProgramOptions(const std::string &name, const std::string &desc, const std::string &syn="-i [FILE] -o [FILE] [OPTIONS]...")
		{
			code_name_=name;
			description_=desc;
			synopsis_=code_name_+" "+syn;
			initialize();
		};

		/**
		 * \brief Destructor of the ProgramOptions object.
		*/
		~ProgramOptions() {};

		/**
		 * \brief Initialize some attribute of the ProgramOptions object.
		*/
		void initialize()
		{
			comment_sign={"#","$"};
			types[""]="";
			types[typeid(bool).name()]="boolean";
			types[typeid(int).name()]="integer";
			types[typeid(size_t).name()]="unsigned integer";
			types[typeid(float).name()]="float";
			types[typeid(double).name()]="double";
			types[typeid(long double).name()]="long double";
			types[typeid(std::string).name()]="string";
			par[MultiKey(std::vector<std::string>({"-parameters","-p"}))]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > ("Parameter file path.", std::vector<std::string>({typeid(std::string).name()}), std::vector<void *>({&param_file}), 1, false, false);

		};

		/**
		 * \brief Formats the displayed informations about parameters.
		 * \param s The string to format.
		 * \param pad Number of tab paddings.
		 * \param N Number of characters in the column.
		 * \return a formated string of characters.
		*/
		std::string format_description(const std::string &s, const size_t &pad=0, const size_t &N=70) const
		{
			std::string tmp(s);
			for (size_t i = 0; i < pad; ++i)
			{
				tmp.insert(0, 1, '\t');
			}
			if (tmp.size()>N)
			{
				size_t pos=N;
				while (pos<tmp.size())
				{
					pos=tmp.find_last_of(" ", pos);
					if (pad>0)
					{
						tmp.replace(pos, 1, 1, '\t');
					}
					for (size_t i = 1; i < pad; ++i)
					{
						tmp.insert(pos, 1, '\t');
					}
					tmp.insert(pos, 1, '\n');
					pos+=N;
				}
			}
			return tmp;
		}

		/**
		 * \brief Checks if an argument is already present in the parameter list.
		 * \param arg The string to check.
		 * \return True or False.
		*/
		bool already_created(const std::string &arg)
		{
			if ( find(MultiKey(arg)) != par_end() )
			{
				#ifdef DEBUG
				std::cout<<"WARNING : The parameter \""<<arg<<"\" was already created."<<std::endl;
				#endif
				return true;
			}
			else
				return false;
		}

		/**
		 * \brief Checks if arguments are already present in the parameter list.
		 * \param arg A vector of strings to check.
		 * \return True or False.
		*/
		bool already_created(std::vector<std::string> &arg)
		{
			if ( find(MultiKey(arg)) != par_end() )
			{
				#ifdef DEBUG
				std::cout<<"WARNING : The parameter \""<<arg[0]<<"\" was already created."<<std::endl;
				#endif
				return true;
			}
			else
				return false;
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		void add_param(const std::string &arg, const std::string &desc, std::vector<std::string> vec_types, std::vector<void *> vec_ref, size_t N=size_t(-1), bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, vec_types, vec_ref, std::min(N,vec_ref.size()), needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		void add_void_param(std::string arg, std::string desc, bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>(), std::vector<void *>(), 0, needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class TT>
		void add_param(const std::string &arg, const std::string &desc, TT &a, size_t N=size_t(-1), bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name()}), std::vector<void *>({&a}), std::min(N,size_t(1)), needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		void add_param(const std::string &arg, const std::string &desc, bool &a, size_t N=0, bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name()}), std::vector<void *>({&a}), std::min(N,size_t(1)), needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class TT>
		void add_param(const std::string &arg, const std::string &desc, std::vector<TT> &a, size_t N=1, bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({"vector",typeid(TT).name()}), std::vector<void *>({&a}), N, needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, size_t N=size_t(-1), bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name()}), std::vector<void *>({&a, &b}), std::min(N,size_t(2)), needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, std::vector<Tb> &b, size_t N=1, bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), "vector", typeid(Tb).name()}), std::vector<void *>({&a, &b}), N, needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, size_t N=size_t(-1), bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name()}), std::vector<void *>({&a, &b, &c}), std::min(N,size_t(3)), needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, std::vector<Tc> &c, size_t N=1, bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), "vector", typeid(Tc).name()}), std::vector<void *>({&a, &b, &c}), N, needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, Td &d, size_t N=size_t(-1), bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), typeid(d).name()}), std::vector<void *>({&a, &b, &c, &d}), std::min(N,size_t(4)), needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, std::vector<Td> &d, size_t N=1, bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), "vector", typeid(d).name()}), std::vector<void *>({&a, &b, &c, &d}), N, needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td, class Te>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, Td &d, Te &e, size_t N=size_t(-1), bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), typeid(d).name(), typeid(e).name()}), std::vector<void *>({&a, &b, &c, &d, &e}), std::min(N,size_t(5)), needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td, class Te>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, Td &d, std::vector<Te> &e, size_t N=1, bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[MultiKey(arg)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), typeid(d).name(), "vector", typeid(e).name()}), std::vector<void *>({&a, &b, &c, &d, &e}), N, needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td, class Te, class Tf>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, Td &d, Te &e, Tf &f, size_t N=size_t(-1), bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[arg]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), typeid(d).name(), typeid(e).name(), typeid(f).name()}), std::vector<void *>({&a, &b, &c, &d, &e, &f}), std::min(N,size_t(6)), needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td, class Te, class Tf>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, Td &d, Te &e, std::vector<Tf> &f, size_t N=1, bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[arg]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), typeid(d).name(), typeid(e).name(), "vector", typeid(f).name()}), std::vector<void *>({&a, &b, &c, &d, &e, &f}), N, needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td, class Te, class Tf, class Tg>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, Td &d, Te &e, Tf &f, Tg &g, size_t N=size_t(-1), bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[arg]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), typeid(d).name(), typeid(e).name(), typeid(f).name(), typeid(g).name()}), std::vector<void *>({&a, &b, &c, &d, &e, &f, &g}), std::min(N,size_t(7)), needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td, class Te, class Tf, class Tg>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, Td &d, Te &e, Tf &f, std::vector<Tg> &g, size_t N=1, bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[arg]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), typeid(d).name(), typeid(e).name(), typeid(f).name(), "vector", typeid(g).name()}), std::vector<void *>({&a, &b, &c, &d, &e, &f, &g}), N, needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td, class Te, class Tf, class Tg, class Th>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, Td &d, Te &e, Tf &f, Tg &g, Th &h, size_t N=size_t(-1), bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[arg]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), typeid(d).name(), typeid(e).name(), typeid(f).name(), typeid(g).name(), typeid(h).name()}), std::vector<void *>({&a, &b, &c, &d, &e, &f, &g, &h}), std::min(N,size_t(8)), needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td, class Te, class Tf, class Tg, class Th>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, Td &d, Te &e, Tf &f, Tg &g, std::vector<Th> &h, size_t N=1, bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[arg]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), typeid(d).name(), typeid(e).name(), typeid(f).name(), typeid(g).name(), "vector", typeid(h).name()}), std::vector<void *>({&a, &b, &c, &d, &e, &f, &g, &h}), N, needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td, class Te, class Tf, class Tg, class Th, class Ti>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, Td &d, Te &e, Tf &f, Tg &g, Th &h, Ti &i, size_t N=size_t(-1), bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[arg]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), typeid(d).name(), typeid(e).name(), typeid(f).name(), typeid(g).name(), typeid(h).name(), typeid(i).name()}), std::vector<void *>({&a, &b, &c, &d, &e, &f, &g, &h, &i}), std::min(N,size_t(9)), needed, false);
			}
		}

		/**
		 * \brief Functions to add parameters to the map.
		*/
		template <class Ta, class Tb, class Tc, class Td, class Te, class Tf, class Tg, class Th, class Ti>
		void add_param(const std::string &arg, const std::string &desc, Ta &a, Tb &b, Tc &c, Td &d, Te &e, Tf &f, Tg &g, Th &h, std::vector<Ti> &i, size_t N=1, bool needed=false)
		{
			if ( !already_created(arg) )
			{
				par[arg]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > (desc, std::vector<std::string>({typeid(a).name(), typeid(b).name(), typeid(c).name(), typeid(d).name(), typeid(e).name(), typeid(f).name(), typeid(g).name(), typeid(h).name(), "vector", typeid(i).name()}), std::vector<void *>({&a, &b, &c, &d, &e, &f, &g, &h, &i}), N, needed, false);
			}
		}

		/**
		 * \brief Copy an existing parameter to another one with different key.
		 * \param arg The argument string of the new argument.
		 * \param existing_arg The argument string of the already defined argument that must be copied.
		*/
		void copy_param(const std::string &arg, const std::string &existing_arg)
		{
			if ( find(existing_arg) != par_end() )
			{
				std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > tmp(find(existing_arg)->second);
				MultiKey a(std::vector<std::string>(((find(existing_arg))->first).keys));
				a.push_back(arg);
				par.erase(existing_arg);
				par[a]=tmp;
			}
			#ifdef DEBUG
			else
				std::cout<<"WARNING : The parameter \""<<existing_arg<<"\" could not be copied because it was not found in the parameter list."<<std::endl;
			#endif
		}

		/**
		 * \brief Copy an existing parameter to another one with different key.
		 * \param arg The argument string of the new argument.
		 * \param existing_arg The argument string of the already defined argument that must be copied.
		 * \param N Minimal number of arguments that must be passed to this parameter.
		*/
		void copy_param(const std::string &arg, const std::string &existing_arg, const size_t &N)
		{
			if ( find(existing_arg) != par_end() )
			{
				std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > tmp(find(existing_arg)->second);
				MultiKey a(std::vector<std::string>(((find(existing_arg))->first).keys));
				a.push_back(arg);
				std::get<3>(tmp)=std::min(N,std::get<2>(tmp).size());
				par.erase(existing_arg);
				par[a]=tmp;
			}
			#ifdef DEBUG
			else
				std::cout<<"WARNING : The parameter \""<<existing_arg<<"\" could not be copied because it was not found in the parameter list."<<std::endl;
			#endif
		}

		/**
		 * \brief Create a new parameter that disables another one: creates "-no_p" from "-p".
		 * \param existing_arg The argument string of the already defined argument that must be copied.
		*/
		void anti_param(const std::string &existing_arg)
		{
			if ( find(existing_arg) != par_end() )
			{
				std::vector<std::string> v((find(existing_arg)->first).keys);
				for (size_t i = 0; i < v.size(); ++i)
				{
					v[i].insert(1,"no_");
				}
				if ( find(MultiKey(v)) != par_end() )
				{
					#ifdef DEBUG
					std::cout<<"WARNING : The anti-parameter of \""<<existing_arg<<"\" was already created."<<std::endl;
					#endif
				}
				else
					par[MultiKey(v)]=std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > ("Cancel the \"" + existing_arg + "\" parameter.", std::get<1>(find(existing_arg)->second), std::get<2>(find(existing_arg)->second), 0, false, false);
			}
			#ifdef DEBUG
			else
				std::cout<<"WARNING : The anti-parameter of \""<<existing_arg<<"\" could not be created because this was not found in the parameter list."<<std::endl;
			#endif
		}

		/**
		 * \brief Displays a man page.
		*/
		void usage()
		{
			std::map<bool,std::string> needed_optional;
			needed_optional[false]="Optional";
			needed_optional[true]="Needed";

			std::cout << std::endl;
			std::cout << "========== USAGE ==========" << std::endl;
			std::cout << std::endl;
			std::cout << std::endl;
			std::cout << "NAME :" << std::endl;
			std::cout << this->format_description(this->code_name_ + " - " + this->description_, 1) << std::endl;
			std::cout << std::endl;
			std::cout << "SYNOPSIS :" << std::endl;
			std::cout << this->format_description(this->synopsis_, 1) << std::endl;
			std::cout << std::endl;
			std::cout << "PARAMETERS :" << std::endl;
			std::cout << this->format_description("Description of all parameters in the following format : ", 1) << std::endl << std::endl;
			std::cout << this->format_description("-parameter <types of argument 1, type of argument 2, ... (types of non needed argument 1, type of non needed argument 2, ...) > : ", 1) << std::endl;
			std::cout << this->format_description("Parameter description.", 2) << std::endl << std::endl;
			std::cout << this->format_description("The parameters are compared case-insensitively and only from the beginning to the first space found. It is thus possible to comment an argument insid a parameter file.", 1) << std::endl;
			std::cout << this->format_description("The types of non needed arguments are in parenthesis.", 1) << std::endl;
			std::cout << this->format_description("The parameter file should have the following structure :", 1) << std::endl;
			std::cout << std::endl;
			std::cout << this->format_description("-parameter1 <- all caracters after that space are not considered", 2) << std::endl;
			std::cout << this->format_description("value1", 2) << std::endl;
			std::cout << this->format_description("$-commented_parameter1", 2) << std::endl;
			std::cout << this->format_description("value1", 2) << std::endl;
			std::cout << this->format_description("#-commented_parameter2", 2) << std::endl;
			std::cout << this->format_description("#commented_value1", 2) << std::endl;
			std::cout << this->format_description("value2", 2) << std::endl;
			std::cout << std::endl;
			std::cout << this->format_description("The \"-parameters\", \"-p\" parameters are reserved for the parameter file path.", 1) << std::endl;
			std::cout << this->format_description("The \"-h\", \"-help\", \"--help\", parameters display user manual.", 1) << std::endl;
			std::cout << std::endl;
			for( std::map<MultiKey, std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > >::const_iterator ii=this->par.begin(); ii!=this->par.end(); ++ii)
			{
				if ( ! std::get<1>(ii->second).empty() )
				{
					std::cout << std::endl;
					std::cout << this->format_description((ii->first).front(), 1) << std::flush;
					if ( std::get<3>(ii->second)!=0 )
					{
						std::cout << " <" << std::flush;
						bool with_bool=false;
						for (size_t i = 0; i < std::get<1>(ii->second).size(); ++i)
						{
							if ( i==0 && ! std::get<1>(ii->second).front().compare(typeid(bool).name()) )
							{
								with_bool=true;
								continue;
							}
							if ( ( i>0 && !with_bool ) || ( i>1 && with_bool ) )
							{
								if ( i>with_bool && i!=std::get<2>(ii->second).size() )
								{
									std::cout << ", " << std::flush;
								}
							}
							if ( i==std::get<3>(ii->second) )
							{
								std::cout << "(" << std::flush;
							}
							if ( !std::get<1>(ii->second).at(i).compare("vector") )
							{
								std::cout << "vector of " << this->types.at(std::get<1>(ii->second).at(i+1)) << std::flush;
								break;
							}
							else
							{
								std::cout << this->types.at(std::get<1>(ii->second).at(i)) << std::flush;
							}
						}
						if ( std::get<2>(ii->second).size()!=std::get<3>(ii->second) )
							std::cout << ")" << std::flush;
						std::cout << ">" << std::flush;
					}
					std::cout << " : (" << needed_optional[std::get<4>(ii->second)] << ")" << std::flush;
					std::cout << std::endl;
				}
				else
				{
					std::cout << this->format_description((ii->first)[0] + " :", 1) << std::endl;
				}
				for (size_t i = 1; i < ii->first.size(); ++i)
				{
					std::cout << this->format_description((ii->first)[i] + " : (Synonym)", 1) << std::endl;
				}
				std::cout << this->format_description(std::get<0>(ii->second), 2) << std::endl;
			}
			std::cout << std::endl;
			std::cout << "AUTHOR :" << std::endl;
			std::cout << this->format_description(author_, 1) << std::endl;
			std::cout << std::endl << std::endl;
		}

		/**
		 * \brief Error displayed if the read argument is not valid.
		 * \param argument The parameter string.
		 * \param n Number of arguments that should have been passed to the parameter.
		 * \param disp Trigger to display or not the error.
		 * \return True or False.
		*/
		bool error_param(const std::string &argument, const size_t n=1, const bool disp=true) const
		{
			if (disp)
			{
				std::cout << std::endl;
				std::cout << "WARNING !!!" << std::endl;
				std::cout << "\""<< argument << "\" parameter needs " << n << " argument(s)." << std::endl;
				std::cout << std::endl;
			}
			return false;
		}

		/**
		 * \brief Checks if a string is commented.
		 * \param str The string to check.
		 * \return True or False.
		*/
		bool check_comment(const std::string &str) const
		{
			bool comment=false;
			for (size_t i = 0; i < comment_sign.size(); ++i)
			{
				if (str.compare(0,comment_sign.at(i).size(),comment_sign.at(i))==0)
				{
					comment=true;
					break;
				}
			}
			return comment;
		}

		/**
		 * \brief Checks if some strings are valid parameter.
		 * \param vec A vector of strings to check.
		 * \param i Index of the first element to check.
		 * \param n Number of elements to check.
		 * \param disp Displays a message if an error is detected.
		 * \return True or False.
		*/
		bool check_param(const std::vector<std::string> &vec, const size_t &i, const size_t n=1, const bool disp=true)
		{
			if (vec.size()<=i+n)
				return error_param(vec.at(i),n,disp);
			for (size_t k = 1; k <= n; ++k)
			{
				if ( check_comment(vec.at(i+k)) || (vec.at(i+k).compare(0,1,"-")==0 && isalpha(vec.at(i+k)[1])!=0) )
					return error_param(vec.at(i),n,disp);
			}
			return true;
		}

		/**
		 * \brief Converts a vector of strings to parameters.
		 * \param vec_ A vector of strings to check.
		 * \param check Triggers the validation of parameters.
		*/
		void param_to_vars(const std::vector<std::string> &vec_, bool check=false)
		{
			std::map<std::string,bool> true_false;
			true_false["1"]=true;
			true_false["0"]=false;
			true_false["yes"]=true;
			true_false["no"]=false;
			true_false["true"]=true;
			true_false["false"]=false;

			std::vector<std::string> vec;
			for (size_t i = 0; i < vec_.size(); ++i)
			{
				if ( !check_comment(vec_.at(i)) )
					vec.push_back(vec_.at(i));
			}
			for (size_t i = 0; i < vec.size(); ++i)
			{
				par_iterator ii=find(vec.at(i));
				if ( ii != par_end() && std::get<5>(ii->second) )
				{
					std::cout << std::endl;
					std::cout << "WARNING !!!" << std::endl;
					std::cout << "The parameter \""  << ii->first.front() << "\" was given several times. Only the first one is considered." << std::endl;
					if ( ii->first.size()>1 )
					{
						std::cout << "Note that this parameter has these synonyms that are considered as equal :" << std::endl;
						for (size_t i = 1; i < ii->first.size(); ++i)
						{
							std::cout << this->format_description((ii->first)[i], 1) << std::endl;
						}
					}
					std::cout << "The current values are : "  << std::endl;
					for (size_t n = 0; n < std::get<2>(ii->second).size(); ++n)
					{
						if ( ! std::get<1>(ii->second)[n].compare(typeid(bool).name()) )
						{
							std::map<bool, std::string> map_bool;
							map_bool[true]="True";
							map_bool[false]="False";
							std::cout << "\t" << map_bool.at(*static_cast<bool*>(std::get<2>(ii->second)[n])) << std::endl;
						}

						else if ( ! std::get<1>(ii->second)[n].compare(typeid(int).name()) )
						{
							std::cout << "\t" << *static_cast<int*>(std::get<2>(ii->second)[n]) << std::endl;
						}

						else if ( ! std::get<1>(ii->second)[n].compare(typeid(size_t).name()) )
						{
							std::cout << "\t" << *static_cast<size_t*>(std::get<2>(ii->second)[n]) << std::endl;
						}

						else if ( ! std::get<1>(ii->second)[n].compare(typeid(float).name()) )
						{
							std::cout << "\t" << *static_cast<float*>(std::get<2>(ii->second)[n]) << std::endl;
						}

						else if ( ! std::get<1>(ii->second)[n].compare(typeid(double).name()) )
						{
							std::cout << "\t" << *static_cast<double*>(std::get<2>(ii->second)[n]) << std::endl;
						}

						else if ( ! std::get<1>(ii->second)[n].compare(typeid(long double).name()) )
						{
							std::cout << "\t" << *static_cast<long double*>(std::get<2>(ii->second)[n]) << std::endl;
						}

						else if ( ! std::get<1>(ii->second)[n].compare(typeid(std::string).name()) )
						{
							std::cout << "\t" << *static_cast<std::string*>(std::get<2>(ii->second)[n]) << std::endl;
						}
					}
					std::cout << std::endl;
					continue;
				}
				if ( ii != par_end() && !std::get<5>(ii->second) )
				{
					size_t nb=0;
					size_t nb_par=std::get<3>(ii->second);
					if ( nb_par>0 && ! std::get<1>(ii->second).front().compare(typeid(bool).name()) )
						--nb_par;
					if ( check_param(vec, i, nb_par) )
					{
						nb=nb_par;
						while ( check_param(vec, i , nb+1, false) )
							++nb;
					}
					size_t j = 0;
					if ( ! std::get<1>(ii->second).front().compare(typeid(bool).name()) )
					{
						++j;
						++nb;
						if ( ii->first.compare(0,4,"-no_")!=0 )
							for (size_t jj = 0; jj < std::get<2>(ii->second).size(); ++jj)
							{
								if ( !std::get<1>(ii->second)[jj].compare(typeid(bool).name()) )
									*static_cast<bool*>(std::get<2>(ii->second)[jj])=true;
								std::get<5>(ii->second)=true;
							}
						else
							for (size_t jj = 0; jj < std::get<2>(ii->second).size(); ++jj)
							{
								if ( !std::get<1>(ii->second)[jj].compare(typeid(bool).name()) )
									*static_cast<bool*>(std::get<2>(ii->second)[jj])=false;
								std::get<5>(ii->second)=true;
							}
					}
					for (size_t k=0; j < nb; ++j,++k)
					{
						if ( ! std::get<1>(ii->second).at(std::max(0,(int)std::get<1>(ii->second).size()-2)).compare("vector") && j>=std::get<1>(ii->second).size()-2)
						{
							size_t ind=std::get<1>(ii->second).size()-2;

							if ( ! std::get<1>(ii->second)[ind+1].compare(typeid(bool).name()) )
							{
								if ( j==0 && !std::get<5>(ii->second) )
									(static_cast< std::vector < bool > * >(std::get<2>(ii->second).back()))->clear();
								(static_cast< std::vector < bool > * >(std::get<2>(ii->second).back()))->push_back(true_false.at(vec[i+k+1].c_str()));
							}

							else if ( ! std::get<1>(ii->second)[ind+1].compare(typeid(int).name()) )
							{
								if ( j==0 && !std::get<5>(ii->second) )
									(static_cast< std::vector < int > * >(std::get<2>(ii->second).back()))->clear();
								(static_cast< std::vector < int > * >(std::get<2>(ii->second).back()))->push_back(atoi(vec[i+k+1].c_str()));
							}

							else if ( ! std::get<1>(ii->second)[ind+1].compare(typeid(size_t).name()) )
							{
								if ( j==0 && !std::get<5>(ii->second) )
									(static_cast< std::vector < size_t > * >(std::get<2>(ii->second).back()))->clear();
								(static_cast< std::vector < size_t > * >(std::get<2>(ii->second).back()))->push_back(static_cast<size_t>(atoi(vec[i+k+1].c_str())));
							}

							else if ( ! std::get<1>(ii->second)[ind+1].compare(typeid(float).name()) )
							{
								if ( j==0 && !std::get<5>(ii->second) )
									(static_cast< std::vector < float > * >(std::get<2>(ii->second).back()))->clear();
								(static_cast< std::vector < float > * >(std::get<2>(ii->second).back()))->push_back(static_cast<float>(atof(vec[i+k+1].c_str())));
							}

							else if ( ! std::get<1>(ii->second)[ind+1].compare(typeid(double).name()) )
							{
								if ( j==0 && !std::get<5>(ii->second) )
									(static_cast< std::vector < double > * >(std::get<2>(ii->second).back()))->clear();
								(static_cast< std::vector < double > * >(std::get<2>(ii->second).back()))->push_back(static_cast<double>(atof(vec[i+k+1].c_str())));
							}

							else if ( ! std::get<1>(ii->second)[ind+1].compare(typeid(long double).name()) )
							{
								if ( j==0 && !std::get<5>(ii->second) )
									(static_cast< std::vector < long double > * >(std::get<2>(ii->second).back()))->clear();
								(static_cast< std::vector < long double > * >(std::get<2>(ii->second).back()))->push_back(static_cast<long double>(atof(vec[i+k+1].c_str())));
							}

							else if ( ! std::get<1>(ii->second)[ind+1].compare(typeid(std::string).name()) )
							{
								if ( j==0 && !std::get<5>(ii->second) )
									(static_cast< std::vector < std::string > * >(std::get<2>(ii->second).back()))->clear();
								(static_cast< std::vector < std::string > * >(std::get<2>(ii->second).back()))->push_back(vec[i+k+1]);
							}
						}
						else
						{
							if ( j>=std::get<2>(ii->second).size() )
								break;
							if ( std::get<1>(ii->second)[j].empty() || ! std::get<1>(ii->second)[j].compare(typeid(bool).name()) )
							{
								try
								{
									*static_cast<bool*>(std::get<2>(ii->second)[j])=true_false.at(vec[i+k+1].c_str());
								}
								catch(const std::out_of_range &oor)
								{
									std::cout << std::endl;
									std::cout << "WARNING !!!" << std::endl;
									std::cout << "The argument of a boolean parameter must be \"0\", \"no\", \"false\", \"1\", \"yes\" or \"true\"." << std::endl;
									std::cout << "The parameter \""  << ii->first.front() << "\" recieved \"" << vec[i+k+1] << "\" as argument." << std::endl;
									if ( ii->first.size()>1 )
									{
										std::cout << "Note that this parameter has synonyms that are considered as equal :" << std::endl;
										for (size_t i = 1; i < ii->first.size(); ++i)
										{
											std::cout << this->format_description((ii->first)[i] + " :", 1) << std::endl;
										}
									}
									std::cout << std::endl;
									continue;
								}
								if ( ii->first.compare(0,4,"-no_")==0 )
									*static_cast<bool*>(std::get<2>(ii->second)[j])=!*static_cast<bool*>(std::get<2>(ii->second)[j]);
							}

							else if ( ! std::get<1>(ii->second)[j].compare(typeid(int).name()) )
							{
								*static_cast<int*>(std::get<2>(ii->second)[j])=atoi(vec[i+k+1].c_str());
							}

							else if ( ! std::get<1>(ii->second)[j].compare(typeid(size_t).name()) )
							{
								*static_cast<size_t*>(std::get<2>(ii->second)[j])=static_cast<size_t>(atoi(vec[i+k+1].c_str()));
							}

							else if ( ! std::get<1>(ii->second)[j].compare(typeid(float).name()) )
							{
								*static_cast<float*>(std::get<2>(ii->second)[j])=static_cast<float>(atof(vec[i+k+1].c_str()));
							}

							else if ( ! std::get<1>(ii->second)[j].compare(typeid(double).name()) )
							{
								*static_cast<double*>(std::get<2>(ii->second)[j])=atof(vec[i+k+1].c_str());
							}

							else if ( ! std::get<1>(ii->second)[j].compare(typeid(long double).name()) )
							{
								*static_cast<long double*>(std::get<2>(ii->second)[j])=atof(vec[i+k+1].c_str());
							}

							else if ( ! std::get<1>(ii->second)[j].compare(typeid(std::string).name()) )
							{
								*static_cast<std::string*>(std::get<2>(ii->second)[j])=vec[i+k+1];
							}

							// else if ( ! std::get<1>(ii->second)[j].compare("function") )
							// {
							// 	*static_cast<std::string*>(std::get<2>(ii->second)[j])=vec[i+k+1];
							// }
						}
						std::get<5>(ii->second)=true;
					}
				}
				else
				{
					if ( check_comment(vec.at(i)) || (vec.at(i).compare(0,1,"-")==0 && isalpha(vec.at(i)[1])!=0) )
					{
						std::cout << std::endl;
						std::cout << "WARNING !!!" << std::endl;
						std::cout << "\"" << vec.at(i).substr(0,vec.at(i).find_first_of(" ")) << "\"" << " is not a valid parameter." << std::endl;
						std::cout << std::endl;
					}
				}
			}
			if (check)
			{
				check_done();
			}
		}

		/**
		 * \brief Reads the parameters from the command line.
		 * \param argc Number of strings.
		 * \param argv A vector of strings.
		 * \param check Triggers the validation of parameters.
		*/
		void param_to_vars(const int argc, const char *argv[], bool check=true)
		{
			this->param_to_vars( std::vector<std::string>(argv+1, argv + argc ), check );
		}

		/**
		 * \brief Checks the parameters stored.
		 * \return True of False.
		*/
		bool check_done()
		{
			bool ok=true;
			std::vector<std::string> invalid;
			for( std::map<MultiKey, std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > >::const_iterator ii=this->par.begin(); ii!=this->par.end(); ++ii)
			{
				if ( std::get<4>(ii->second) && ! std::get<5>(ii->second) )
				{
					#ifdef DEBUG
					std::cout<<"Invalid : "<<ii->first.front()<<" : "<<std::get<4>(ii->second)<<" "<<std::get<5>(ii->second)<<std::endl;
					if ( ii->first.size()>1 )
					{
						std::cout << "Note that this parameter has synonyms that are considered as equal :" << std::endl;
						for (size_t i = 1; i < ii->first.size(); ++i)
						{
							std::cout << this->format_description((ii->first)[i] + " :", 1) << std::endl;
						}
					}
					#endif
					invalid.push_back(ii->first.front());
					ok=false;
				}
			}
			if (!ok)
			{
				std::cout<<std::endl;
				std::cout<<"******************************************"<<std::endl;
				std::cout<<"******* /!\\ Invalid Parameters /!\\********"<<std::endl;
				std::cout<<"******************************************"<<std::endl;
				for (size_t i = 0; i < invalid.size(); ++i)
				{
					std::cout<<invalid.at(i).substr(0,invalid.at(i).find_first_of(" "))<<" needs at least "<<std::get<3>(par[invalid.at(i)])<<" argument(s)."<<std::endl;
				}
				std::cout<<"******************************************"<<std::endl;
				usage();
			}
			return ok;
		}

		/**
		 * \brief Reads parameters from command line and file ("-p" parameter).
		 * \param vec A vector of strings.
		 * \param check Triggers the validation of parameters.
		 * \return True of False.
		*/
		bool read(const std::vector<std::string> &vec, const bool check=true)
		{
			std::string str;
			for (size_t i = 0; i < vec.size(); ++i)
			{
				if ( !vec.at(i).compare("-h") || !vec.at(i).compare("-help") || !vec.at(i).compare("--help"))
				{
					usage();
					exit(0);
				}
			}

			param_to_vars( vec );

			if ( ! param_file.empty() )
			{
				std::vector<std::string> vec_file;
				std::ifstream ifile;
				ifile.open(param_file.c_str());
				if (ifile.fail())
				{
					std::cout<<"***************************************"<<std::endl;
					std::cout<<"Unable to read data file"<<std::endl;
					std::cout<<"Parameters are set to default values"<<std::endl;
					std::cout<<"***************************************"<<std::endl;
					sleep(3);
				}
				while (getline(ifile,str))
				{
					if ( check_comment(str) )
						continue;
					if ( ! (str.compare(0,1,"-")) && isalpha(str[1]) )
						vec_file.push_back(str);
					else
					{
						while (str.size()>0)
						{
							size_t pos=str.find_first_of(" \t");
							vec_file.push_back(str.substr(0,pos));
							if (pos!=std::string::npos)
								str.erase(0, pos+1);
							else
								str.clear();
						}
					}
				}
				ifile.close();
				param_to_vars( vec_file );
			}
			if (check)
			{
				return check_done();
			}
			else
			{
				return true;
			}
		}

		/**
		 * \brief Reads parameters from command line and file ("-p" parameter).
		 * \param argc Number of parameters.
		 * \param argv A char** containing the parameters.
		 * \return True of False.
		*/
		bool read(const int argc, const char *argv[])
		{
			return this->read( std::vector<std::string>(argv+1, argv + argc ) );
		}

		/**
		 * \brief Reads parameters from command line and file ("-p" parameter).
		 * \param argc Number of parameters.
		 * \param argv A char** containing the parameters.
		 * \return True of False.
		*/
		bool read(const int argc, char *argv[])
		{
			return this->read( std::vector<std::string>(argv+1, argv + argc ) );
		}

		/**
		 * \brief Outputs operator to write the given parameters and their arguments.
		 * \param out Incoming stream.
		 * \param a ProgramOptions object to display.
		 * \return A stream.
		*/
		friend std::ostream& operator<<(std::ostream &out, const ProgramOptions &a);

		/**
		 * \brief Writes the parameters to a file.
		 * \param name Path to the file.
		 * \param mode std::ios_base mode.
		*/
		void write_parameters(const std::string name, const std::ios_base::openmode mode=std::ios_base::trunc)
		{
			std::ofstream file;
			file.open(name.c_str(),mode);
			file<<*this;
			file.close();
		};

		/**
		 * \name Iterators
		 * @{
		*/
		typedef std::map<MultiKey, std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > >::iterator par_iterator;
		typedef std::map<MultiKey, std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > >::reverse_iterator par_reverse_iterator;
		typedef std::map<MultiKey, std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > >::const_iterator par_const_iterator;
		typedef std::map<MultiKey, std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > >::const_reverse_iterator par_const_reverse_iterator;

		/**
		 * \brief Returns a read/write iterator to the first parameter.
		*/
		par_iterator par_begin()
		{
			return par.begin();
		}

		/**
		 * \brief Returns a read/write iterator past to the last parameter.
		*/
		par_iterator par_end()
		{
			return par.end();
		}

		/**
		 * \brief Returns a read/write iterator to the last parameter.
		*/
		par_reverse_iterator par_rbegin()
		{
			return par.rbegin();
		}

		/**
		 * \brief Returns a read/write iterator past to the first parameter.
		*/
		par_reverse_iterator par_rend()
		{
			return par.rend();
		}

		/**
		 * \brief Returns a read iterator to the last parameter.
		*/
		par_const_iterator par_cbegin() const
		{
			return par.cbegin();
		}

		/**
		 * \brief Returns a read iterator past to the last parameter.
		*/
		par_const_iterator par_cend() const
		{
			return par.cend();
		}

		/**
		 * \brief Returns a read iterator to the first parameter.
		*/
		par_const_reverse_iterator par_crbegin() const
		{
			return par.crbegin();
		}

		/**
		 * \brief Returns a read iterator past to the first parameter.
		*/
		par_const_reverse_iterator par_crend() const
		{
			return par.crend();
		}
		/**@} End iterators*/

		/**
		 * \brief Checks if a parameter was given or not.
		 * \param a The key to look for.
		 * \return An iterator to the parameter if it was found or map::end() otherwise.
		*/
		par_iterator find(const MultiKey &a)
		{
			par_iterator i=par_begin();
			for (; i != par_end(); ++i)
			{
				if ( a == i->first )
				{
					return i;
				}
			}
			return par_end();
		}

		/**
		 * \brief Checks if a parameter was given or not.
		 * \param a The key to look for.
		 * \return An iterator to the parameter if it was found or map::end() otherwise.
		*/
		par_const_iterator find(const MultiKey &a) const
		{
			par_const_iterator i=par_cbegin();
			for (; i != par_cend(); ++i)
			{
				if ( a == i->first )
				{
					return i;
				}
			}
			return par_cend();
		}

		/**
		 * \brief Checks if a paramters was given or not.
		 * \param a The key to look for.
		 * \return A reverse_iterator to the parameter if it was found or map::rend() otherwise.
		*/
		par_reverse_iterator rfind(const MultiKey &a)
		{
			par_reverse_iterator i=par_rbegin();
			for (; i != par_rend(); ++i)
			{
				if ( a == i->first )
				{
					return i;
				}
			}
			return par_rend();
		}

		/**
		 * \brief Checks if a paramters was given or not.
		 * \param a The key to look for.
		 * \return A const_iterator to the parameter if it was found or map::end() otherwise.
		*/
		par_const_iterator cfind(const MultiKey &a) const
		{
			par_const_iterator i=par_cbegin();
			for (; i != par_cend(); ++i)
			{
				if ( a == i->first )
				{
					return i;
				}
			}
			return par_cend();
		}

		/**
		 * \brief Checks if a paramters was given or not.
		 * \param a The key to look for.
		 * \return A const_iterator to the parameter if it was found or map::rend() otherwise.
		*/
		par_const_reverse_iterator crfind(const MultiKey &a)
		{
			par_const_reverse_iterator i=par_crbegin();
			for (; i != par_crend(); ++i)
			{
				if ( a == i->first )
				{
					return i;
				}
			}
			return par_crend();
		}

		/**
		 * \brief Checks if a paramters is needed or not.
		 * \param a The parameters to check.
		 * \return True or False.
		*/
		bool is_needed(const std::string &a) const
		{
			if ( std::get<4>((find(a)->second)) )
				return true;
			else
				return false;
		}

		/**
		 * \brief Checks if a paramters was given or not.
		 * \param a The key to look for.
		 * \return True or False.
		*/
		bool was_given(const std::string &a) const
		{
			par_const_iterator i=find(a);
			if ( i != this->par.end() )
			{
				if ( std::get<5>(i->second) )
					return true;
				else
					return false;
			}
			else
				return false;
		}

		/**
		 * \brief Gets the name of the code.
		 * \return A reference to the string containing the name of the code.
		*/
		std::string& code_name()
		{
			return this->code_name_;
		}

		/**
		 * \brief Gets the name of the code.
		 * \return A reference to the string containing the name of the code.
		*/
		const std::string& code_name() const
		{
			return this->code_name_;
		}

		/**
		 * \brief Sets the name of the code.
		 * \param a A string containing the name of the code.
		*/
		void code_name(std::string &a)
		{
			this->code_name_=a;
		}

		/**
		 * \brief Gets the synopsis of the code.
		 * \return A reference to the string containing the synopsis of the code.
		*/
		std::string& synopsis()
		{
			return this->synopsis_;
		}

		/**
		 * \brief Gets the synopsis of the code.
		 * \return A reference to the string containing the synopsis of the code.
		*/
		const std::string& synopsis() const
		{
			return this->synopsis_;
		}

		/**
		 * \brief Sets the synopsis of the code.
		 * \param a A string containing the synopsis of the code.
		*/
		void synopsis(std::string &a)
		{
			this->synopsis_=a;
		}

		/**
		 * \brief Gets the description of the code.
		 * \return A reference to the string containing the description of the code.
		*/
		std::string& description()
		{
			return this->description_;
		}

		/**
		 * \brief Gets the description of the code.
		 * \return A reference to the string containing the description of the code.
		*/
		const std::string& description() const
		{
			return this->description_;
		}

		/**
		 * \brief Sets the description of the code.
		 * \param a A string containing the description of the code.
		*/
		void description(std::string &a)
		{
			this->description_=a;
		}

		/**
		 * \brief Gets the parameter file path.
		 * \return A reference to the string containing the parameter file path.
		*/
		std::string& parameter_file_path()
		{
			return this->param_file;
		}

		/**
		 * \brief Gets the parameter file path.
		 * \return A reference to the string containing the parameter file path.
		*/
		const std::string& parameter_file_path() const
		{
			return this->param_file;
		}

		/**
		 * \brief Sets the parameter file path.
		 * \param a A string containing the parameter file path.
		*/
		void parameter_file_path(std::string &a)
		{
			this->param_file=a;
		}

		/**
		 * \brief Gets the author(s).
		 * \return A reference to the string containing the author(s).
		*/
		std::string& author()
		{
			return this->author_;
		}

		/**
		 * \brief Gets the author(s).
		 * \return A reference to the string containing the author(s).
		*/
		const std::string& author() const
		{
			return this->author_;
		}

		/**
		 * \brief Sets the author(s) of the code.
		 * \param a A string containing the author(s) of the code.
		*/
		void author(std::string &a)
		{
			this->author_=a;
		}

	private:
		std::vector<std::string> comment_sign; /**<Strings that are recognized as a comment trigger.*/
		std::string code_name_; /**<String storing the name of the code.*/
		std::string synopsis_; /**<String storing the synopsis of the code.*/
		std::string description_; /**<String storing the description of the code.*/
		std::string author_; /**<String storing the author(s) of the code.*/
		std::string param_file; /**<String storing the parameter file.*/
		std::map<std::string,std::string> types; /**<Map of strings storing the data types.*/
		std::map<MultiKey,std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool> > par; /**<Map storing the parameters and their data.*/
	};

	std::ostream& operator<<(std::ostream &out, const ProgramOptions &a)
	{
		std::map<int,std::string> tf;
		tf[0]="No";
		tf[1]="Yes";
		out<<"********************************************"<<std::endl;
		out<<"************ Given Parameters **************"<<std::endl;
		out<<"********************************************"<<std::endl;
		for (std::map<MultiKey, std::tuple<std::string, std::vector<std::string>, std::vector<void *>, size_t, bool, bool > >::const_iterator i = a.par.begin(); i != a.par.end(); ++i)
		{
			if (std::get<5>((*i).second))
			{
				out<<std::endl;
				for (size_t j = 0; j < (*i).first.keys.size(); ++j)
				{
					out<<(*i).first.keys.at(j);
					if ( (*i).first.keys.size()>1 && j<(*i).first.keys.size()-1 )
					{
						out<<" / ";
					}
				}
				out<<std::endl;

				for (size_t j = 0; j < std::get<2>((*i).second).size(); ++j)
				{
					if ( ! std::get<1>((*i).second).at(std::max(0,(int)std::get<1>((*i).second).size()-2)).compare("vector") && j>=std::get<1>((*i).second).size()-2)
					{

						if ( ! std::get<1>((*i).second).back().compare(typeid(bool).name()) )
						{
							const size_t size=(*static_cast< std::vector < bool > * >(std::get<2>((*i).second).back())).size();
							for (size_t k = 0; k < size; ++k)
								out<<(*static_cast< std::vector < bool > * >(std::get<2>((*i).second).back())).at(k)<<std::endl;
						}

						else if ( ! std::get<1>((*i).second).back().compare(typeid(int).name()) )
						{
							const size_t size=(*static_cast< std::vector < int > * >(std::get<2>((*i).second).back())).size();
							for (size_t k = 0; k < size; ++k)
								out<<(*static_cast< std::vector < int > * >(std::get<2>((*i).second).back())).at(k)<<std::endl;
						}

						else if ( ! std::get<1>((*i).second).back().compare(typeid(size_t).name()) )
						{
							const size_t size=(*static_cast< std::vector < size_t > * >(std::get<2>((*i).second).back())).size();
							for (size_t k = 0; k < size; ++k)
								out<<(*static_cast< std::vector < size_t > * >(std::get<2>((*i).second).back())).at(k)<<std::endl;
						}

						else if ( ! std::get<1>((*i).second).back().compare(typeid(float).name()) )
						{
							const size_t size=(*static_cast< std::vector < float > * >(std::get<2>((*i).second).back())).size();
							for (size_t k = 0; k < size; ++k)
								out<<(*static_cast< std::vector < float > * >(std::get<2>((*i).second).back())).at(k)<<std::endl;
						}

						else if ( ! std::get<1>((*i).second).back().compare(typeid(double).name()) )
						{
							const size_t size=(*static_cast< std::vector < double > * >(std::get<2>((*i).second).back())).size();
							for (size_t k = 0; k < size; ++k)
								out<<(*static_cast< std::vector < double > * >(std::get<2>((*i).second).back())).at(k)<<std::endl;
						}

						else if ( ! std::get<1>((*i).second).back().compare(typeid(long double).name()) )
						{
							const size_t size=(*static_cast< std::vector < long double > * >(std::get<2>((*i).second).back())).size();
							for (size_t k = 0; k < size; ++k)
								out<<(*static_cast< std::vector < long double > * >(std::get<2>((*i).second).back())).at(k)<<std::endl;
						}

						else if ( ! std::get<1>((*i).second).back().compare(typeid(std::string).name()) )
						{
							const size_t size=(*static_cast< std::vector < std::string > * >(std::get<2>((*i).second).back())).size();
							for (size_t k = 0; k < size; ++k)
								out<<(*static_cast< std::vector < std::string > * >(std::get<2>((*i).second).back())).at(k)<<std::endl;
						}

						break;
					}
					else
					{
						if ( std::get<1>((*i).second)[j].empty() || ! std::get<1>((*i).second)[j].compare(typeid(bool).name()) )
						{
							if ( ! std::get<2>((*i).second).empty() )
							{
								if ( (*i).first.compare(0,4,"-no_")==0 )
									out<<tf.at(!*static_cast<bool*>(std::get<2>((*i).second)[j]))<<std::endl;
								else
									out<<tf.at(*static_cast<bool*>(std::get<2>((*i).second)[j]))<<std::endl;
							}
						}

						else if ( ! std::get<1>((*i).second)[j].compare(typeid(int).name()) )
						{
							out<<*static_cast<int*>(std::get<2>((*i).second)[j])<<std::endl;
						}

						else if ( ! std::get<1>((*i).second)[j].compare(typeid(size_t).name()) )
						{
							out<<*static_cast<size_t*>(std::get<2>((*i).second)[j])<<std::endl;
						}

						else if ( ! std::get<1>((*i).second)[j].compare(typeid(float).name()) )
						{
							out<<*static_cast<float*>(std::get<2>((*i).second)[j])<<std::endl;
						}

						else if ( ! std::get<1>((*i).second)[j].compare(typeid(double).name()) )
						{
							out<<*static_cast<double*>(std::get<2>((*i).second)[j])<<std::endl;
						}

						else if ( ! std::get<1>((*i).second)[j].compare(typeid(long double).name()) )
						{
							out<<*static_cast<long double*>(std::get<2>((*i).second)[j])<<std::endl;
						}

						else if ( ! std::get<1>((*i).second)[j].compare(typeid(std::string).name()) )
						{
							out<<*static_cast<std::string*>(std::get<2>((*i).second)[j])<<std::endl;
						}
					}
				}
			}
		}
		out<<std::endl;
		return out;
	}

} // namespace generic

#endif
