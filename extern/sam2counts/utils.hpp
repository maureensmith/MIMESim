#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <string>
#include <cctype> 

namespace utils
{
inline unsigned string_to_unsigned(std::string::iterator& it)
{
    unsigned ret = 0;
    do
    {
        const auto val = *it;
        if (not std::isdigit(val))
        {
            return ret;
        }
        ret *= 10;
        ret += val - '0';
        ++it;
    }
    while (true);
}

template <int N>
std::string::iterator find_tab_in_string(std::string::iterator it)
{
    int i = 0;
    while (i != N)
    {
        if (*it == '\t')
        {
            ++i;
        }
        ++it;
    }
    return it;
}

unsigned long choose(const int n, const int k);

}



#endif

