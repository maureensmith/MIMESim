#ifndef REF_MAP_HPP
#define REF_MAP_HPP

#include <utility>
#include "nucleobase.hpp"

#include <iterator>
#include <cstring>


namespace ref
{
class ref_map final
{
  public:
    using data_type = std::pair<unsigned, nucleotid::nucleobase>;
    using size_type = std::size_t;

    ref_map() = default;

    ~ref_map()
    {
        delete[] data;
    }

    auto begin()
    {
        return data;
    }

    auto begin() const
    {
        return data;
    }

    auto end() const
    {
        return fake_end;
    }

    auto end()
    {
        return fake_end;
    }

    void remove(data_type* pos)
    {
        --fake_end;
        *pos = *fake_end;
    }

    void reserve(const size_type size)
    {
        if (size >= capacity)
        {
            capacity = size + 1;
            auto tmp = data;
            const auto ns = this->size();
            data = new data_type[capacity];
            //TODO: work around: seit gcc8 compile fehler
            std::memcpy(reinterpret_cast<char*> (data), reinterpret_cast<char*> (tmp), sizeof (data_type) * ns);
            fake_end = data + ns;
            delete[] tmp;
        }
    }

    void clear()
    {
        fake_end = data;
    }

    void add(std::pair<unsigned, nucleotid::nucleobase>&& pair)
    //void add(const unsigned p, const char n)
    {
        *fake_end = pair;
        ++fake_end;
        //*(fake_end - 1) = std::make_pair(p, nucleobase{n});
    }

    size_type size() const noexcept
    {
        return std::distance(begin(), end());
    }

  private:
    size_type capacity{0};
    data_type* data{nullptr};
    data_type* fake_end{nullptr};
};
}

#endif

