
// Created by msmith on 15.01.20.
//

#ifndef MIMESIM_COUNTER_H
#define MIMESIM_COUNTER_H

#include <vector>
#include <string>

namespace count {
    using count_type = unsigned;

    class counter_1 final {
    public:

        // Delete default constructor
        counter_1() = delete;

        ~counter_1() = default;


        counter_1(const unsigned int L, const unsigned int q)
                : data(L,std::vector<count_type>(q)) {
        }

        counter_1(const unsigned int L, const unsigned int q, const int initialValue)
                : data(L,std::vector<count_type>(q, initialValue)) {
        }

        void increment(const unsigned pos, const unsigned symbol);

        void decrement(const unsigned pos, const unsigned symbol);

        void count(const unsigned pos, const unsigned symbol, const unsigned times);

        const count_type get(const unsigned pos, const unsigned symbol);

        void write_to_file(const std::string &out_file);

    private:
        //containing counts for all q nucleotides for each position (counted from 0, for input/output -1)
        std::vector <std::vector<count_type>> data;
    };

    class counter_2 final {
    public:
        // Delete default constructor
        counter_2() = delete;

        ~counter_2() = default;


        counter_2 (const unsigned int length, const unsigned int symbols)
        : L(length), q(symbols), data(length*(length-1)/2,std::vector<count_type>(symbols*symbols)) {}

        counter_2 (const unsigned int length, const unsigned int symbols, const int initialValue)
                : L(length), q(symbols), data(length*(length-1)/2,std::vector<count_type>(symbols*symbols, initialValue)) {}

        void increment(const unsigned pos1, const unsigned pos2, const unsigned symbol1, const unsigned symbol2);

        void decrement(const unsigned pos1, const unsigned pos2, const unsigned symbol1, const unsigned symbol2);

        void count(const unsigned pos1, const unsigned pos2, const unsigned symbol1, const unsigned symbol2, const int times);

        /*
         * Fill value for all positions pairs with pos containing symbol and the second positions being wildtype =1
         */
        void count(const unsigned pos, const unsigned symbol, const int times);


        void write_to_file(const std::string &out_file);

    private:
        unsigned L;
        unsigned q;
        //containing counts for all q nucleotides for each position (counted from 0, for input/output -1)
        std::vector <std::vector<count_type>> data;
        unsigned getPositionIndex(const unsigned pos1, const unsigned pos2);
        unsigned getSymbolIndex(const unsigned symbol1, const unsigned symbol2);
    };


    /**
     * Contains 1D and 2D counter for bound and unbound fractions
     */
    class counter_collection final {
    public:
        // Delete default constructor
        counter_collection() = delete;

        ~counter_collection() = default;

        counter_collection (const unsigned int L, const unsigned int q)
                : counter_bound_1d(L, q), counter_unbound_1d(L, q), counter_bound_2d(L, q), counter_unbound_2d(L, q) {}

        counter_1 counter_bound_1d;
        counter_1 counter_unbound_1d;

        counter_2 counter_bound_2d;
        counter_2 counter_unbound_2d;

    };
}

#endif //MIMESIM_COUNTER_H
