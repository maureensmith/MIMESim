#ifndef NUCLEOBASE_HPP_
#define NUCLEOBASE_HPP_

#include <string>
#include <array>
#include <cstdint>


namespace nucleotid
{
class nucleobase final
{
  private:
    static const std::array<char, 6> bases;
    uint8_t base{5};
    //uint8_t bla[11] = {0, 1, 0, 2, 0, 0, 0, 0, 0, 3, 3};

  public:
    nucleobase() = default;
    ~nucleobase() = default;

    // standard copy constructor
    nucleobase(const nucleobase&) = default;
    // standard assignment operator
    nucleobase& operator= (const nucleobase&) = default;

    explicit nucleobase(char base);
    explicit nucleobase(const int base_id);

    char get() const noexcept
    {
        return nucleobase::bases[base];
    }

    uint8_t to_id() const
    {
        return base;
    }

    bool operator==(const nucleobase& rhs) const noexcept
    {
        return to_id() == rhs.to_id();
    }
};
}

#endif
