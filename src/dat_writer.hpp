/**
 * \file dat_writer.hpp
 *
 * \brief Provides a class to write DAT files.
 *
*/
#ifndef DAT_WRITER_HPP
#define DAT_WRITER_HPP
#include <fstream>
#include <memory>

namespace generic
{

    /**
     * \ingroup generic
     * \brief This is a class that provides methods to write DAT files.
    */
    class DatWriter
    {
    public:

        /**
         * \brief Constructor of the DatWriter.
         * \param file_name Name of the file.
         * \param mode Open mode. The default mode is std::ofstream::out | std::ofstream::binary | std::ofstream::trunc.
         * \param delimiter Delimiter used to create DAT file columns. The default value is a comma.
        */
        DatWriter(const std::string& file_name, const std::ios_base::openmode& mode = std::ofstream::out | std::ofstream::binary | std::ofstream::trunc, const std::string& delimiter = ",") : _file_name(file_name), _mode(mode), _delimiter(delimiter)
        {
            open();
        };

        /**
         * \brief Constructor of the DatWriter.
         * \param file_name Name of the file.
         * \param delimiter Delimiter used to create DAT file columns. The default value is a comma.
        */
        DatWriter(const std::string& file_name, const std::string& delimiter) : _file_name(file_name), _delimiter(delimiter)
        {
            _mode = std::ofstream::out | std::ofstream::binary | std::ofstream::trunc;
            open();
        };

        /**
         * \brief Constructor of the DatWriter.
        */
        DatWriter()
        {
        };

        /**
         * \brief Constructor of the DatWriter.
        */
        ~DatWriter()
        {
            close();
        };

        /**
         * \brief Provides read/write access to the file name.
         * \return Address of the _file_name attribute.
        */
        std::string& file_name()
        {
            return _file_name;
        }

        /**
         * \brief Provides read/write access to the opening mode.
         * \return Address of the _mode attribute.
        */
        std::ios_base::openmode& mode()
        {
            return _mode;
        }

        /**
         * \brief Provides read/write access to the delimiter.
         * \return Address of the _delimiter attribute.
        */
        std::string& delimiter()
        {
            return _delimiter;
        }

        /**
         * \brief Opens a new file.
         * \param file_name Name of the file.
         * \param mode Open mode.
         * \param delimiter The delimiter used.
        */
        void open(const std::string& file_name, const std::ios_base::openmode& mode = std::ofstream::out | std::ofstream::binary | std::ofstream::trunc, const std::string& delimiter = ",")
        {
            _file_name = file_name;
            _mode = mode;
            _delimiter = delimiter;
            open();
        }

        /**
         * \brief Opens a new file.
         * \param file_name Name of the file.
         * \param delimiter The delimiter used.
        */
        void open(const std::string& file_name, const std::string& delimiter = ",")
        {
            _file_name = file_name;
            _mode = std::ofstream::out | std::ofstream::binary | std::ofstream::trunc;
            _delimiter = delimiter;
            open();
        }

        /**
         * \brief Closes the file.
        */
        void close()
        {
            _file.close();
        }

        /**
         * \brief Writes the header.
         * \param cols Column names.
        */
        template <typename ...ColNames>
        void write_header(const std::string &title, const size_t &I, const size_t &J, ColNames... cols)
        {
            _n_cols = 0;
            _file<<"TITLE = "<<title<<std::endl;
            _file<<"VARIABLES = ";
            write_header_sub(cols...);
            _file<<"ZONE T = \"Values\", I="<<I<<", J="<<J<<std::endl;
        };

        /**
         * \brief Writes the row data.
         * \param cols Column names.
        */
        template <typename ...ColNames>
        inline
        void write_row(ColNames... cols)
        {
            std::string n = std::to_string(sizeof...(cols));
            if (sizeof...(cols) != _n_cols)
            {
                throw std::runtime_error("Error in DatWriter: The number of elements must be equal to number of columns : " + n + " != "+std::to_string(_n_cols));
            }
            else
            {
                write_row_sub(cols...);
            }
        };

    protected:
        /**
         * \brief Opens a new file using the parameters given by the attributes.
        */
        void open()
        {
            if (_file.is_open())
            {
                close();
            }
            _file.open(_file_name, _mode);
        }

        /**
         * \brief Writes the header.
         * \param col Column name currently processed.
         * \param cols Column names.
         * \remark This function should be only called by the function write_header.
        */
        template <typename ...ColNames>
        inline
        void write_header_sub(const std::string& col, ColNames... cols)
        {
            _file<<"\""<<col<<"\""<<_delimiter;
            write_header_sub(cols...);
            ++_n_cols;
        };

        /**
         * \brief Writes the header.
         * \param col Column name currently processed.
         * \remark This function should be only called by the function write_header.
        */
        inline
        void write_header_sub(const std::string& col)
        {
            _file<<"\""<<col<<"\""<<std::endl;
            ++_n_cols;
        };

        /**
         * \brief Writes the row data.
         * \param col Column data currently processed.
         * \param cols Other column names.
         * \remark This function should be only called by the function write_header.
        */
        template <typename T, typename ...ColNames>
        inline
        void write_row_sub(const T& col, ColNames... cols)
        {
            _file<<col<<_delimiter;
            write_row_sub(cols...);
        };

        /**
         * \brief Writes the row data.
         * \param col Column data currently processed.
         * \remark This function should be only called by the function write_header.
        */
        template <typename T>
        inline
        void write_row_sub(const T& col)
        {
            _file<<col<<std::endl;
        };

    private:
        std::ofstream _file; /**<Handle to the file.*/
        std::string _file_name; /**<File name.*/
        std::ios_base::openmode _mode; /**<Opening mode of the file.*/
        std::string _delimiter; /**<Delimiter used.*/
        std::size_t _n_cols = 0; /**<Number of columns counted in the header.*/
    };

} // namespace generic

#endif // DAT_WRITER_HPP
