/**
 * \file logger.hpp
 *
 * \brief Provides a structure to log informations.
 *
*/
#ifndef LOGGER_HPP
#define LOGGER_HPP
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>


namespace generic
{

	/**
	 * \ingroup generic
	 * \brief Initialise a boost logger
	 */
	void init_logger(
	    const std::string &log_level=std::string(),
	    const std::string &log_file=std::string(),
	    const std::string &format=std::string("%TimeStamp% - %Severity%\t- %Message%"))
	{
	    boost::log::add_common_attributes();
	    boost::log::register_simple_formatter_factory< boost::log::trivial::severity_level, char >("Severity");

	    if (!log_file.empty())
	    {
	        boost::log::add_file_log(
	            boost::log::keywords::file_name = log_file,
	            boost::log::keywords::format = format);
	    }
	    boost::log::add_console_log(
	        std::cout,
	        boost::log::keywords::format = format);

	    if (log_level == "trace")
	    {
	        boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::trace);
	    }
	    else if (log_level == "debug")
	    {
	        boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::debug);
	    }
	    else if (log_level == "info")
	    {
	        boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::info);
	    }
	    else if (log_level == "warning")
	    {
	        boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::warning);
	    }
	    else if (log_level == "error")
	    {
	        boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::error);
	    }
	    else if (log_level == "fatal")
	    {
	        boost::log::core::get()->set_filter(boost::log::trivial::severity >= boost::log::trivial::fatal);
	    }
	    else if (log_level == "quiet")
	    {
	        boost::log::core::get()->set_filter(boost::log::trivial::severity > boost::log::trivial::fatal);
	    }
	    else
	    {
	        throw std::invalid_argument("The log level argument must be in the following list: [trace, debug, info, warning, error, fatal, quiet].");
	    }
	}

} // namespace generic

#endif // LOGGER_HPP
