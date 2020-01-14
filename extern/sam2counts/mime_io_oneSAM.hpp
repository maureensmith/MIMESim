#ifndef MIME_IO_HPP_
#define MIME_IO_HPP_

#include <string>

namespace mime_io
{
void analyse_positions(const std::string& ref,
                       const std::string& sam_a,
                       const std::string& out_file,
                       const int dimension);
}

#endif
