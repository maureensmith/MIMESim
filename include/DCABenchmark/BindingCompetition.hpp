//
// Created by Smith, Maureen on 04.06.18.
//

#ifndef DCABENCHMARK_BINDINGCOMPETITION_HPP
#define DCABENCHMARK_BINDINGCOMPETITION_HPP

#include "Species.hpp"
#include "Constants.hpp"
#include "../extern/CppNumericalSolvers/cppoptlib/meta.h"
#include "../extern/CppNumericalSolvers/cppoptlib/boundedproblem.h"
#include <valarray>

//template<typename T>
//TODO: Umbau nach counts
class UnboundProtein : public cppoptlib::BoundedProblem<double> {
//class UnboundProtein : public cppoptlib::BoundedProblem<double> {
private:
    /**** Constants regarding ODE solving *****/
    //the total amount of Protein
    static constexpr double B_TOT = 2.0;
    //TODO: Umbau nach counts
    //static constexpr int B_TOT = 12*pow(10,6);

    //species kds
    const std::valarray<double> kds;
    //species frequencies
    const std::valarray<double> frequencies;
    //TODO: Umbau nach counts
    //species counts
    const std::valarray<unsigned int> counts;

    std::valarray<double> getSpeciesKds(const species::species_map& spec);

    std::valarray<double> getSpeciesFrequencies(const species::species_map& spec);

    //TODO: Umbau nach counts
    std::valarray<unsigned int> getSpeciesCounts(const species::species_map& spec);

public:
    //TODO: Umbau nach counts
    //using Superclass = cppoptlib::BoundedProblem<double>;
    using Superclass = cppoptlib::BoundedProblem<double>;
    using typename Superclass::TVector;
    //UnboundProtein(const species::species_map& spec) : Superclass(1), kds(getSpeciesKds(spec)) , frequencies(getSpeciesFrequencies(spec)) {};
    //TODO: Umbau nach counts
    UnboundProtein(const species::species_map& spec) : Superclass(1), kds(getSpeciesKds(spec)),
                                                       frequencies(getSpeciesFrequencies(spec)), counts(getSpeciesCounts(spec)) {};

    //TODO: Umbau nach counts
    // the objective to be minimised
    double value(const TVector &x) {
    //double value(const TVector &x) {
        //TODO: Umbau nach counts
        return   pow(B_TOT - (x[0]*(frequencies/(kds+double(x[0]))).sum())-x[0], 2);
        //return  pow(B_TOT - (x[0]*(counts/(kds+double(x[0]))).sum())-x[0], 2);
    }

    //TODO: Umbau nach counts
    //double solve(std::valarray<double>& S_bound, std::valarray<double>& S_unbound);
    double solve(std::valarray<unsigned int>& S_bound, std::valarray<unsigned int>& S_unbound);

    static const int getB_tot();

    const std::valarray<double> &getKds() const;

    const std::valarray<double> &getFrequencies() const;

};


#endif //DCABENCHMARK_BINDINGCOMPETITION_HPP
