//
// Created by Smith, Maureen on 31.05.18.
//

#ifndef DCABENCHMARK_FUNCTIONALSEQUENCE_HPP
#define DCABENCHMARK_FUNCTIONALSEQUENCE_HPP

#include <vector>
#include <iostream>

class FunctionalSequence {

private:
    // original Kd for each site
    const std::vector<double> kds;
    // original epistasis values for each position pair
    const std::vector<double> epistasis;

    //singleton class has a private constructor which is called by a public method
    FunctionalSequence():kds(drawKdValues()), epistasis(drawEpistasis()){};
    // draws the Kds value for each L postions
    std::vector<double> drawKdValues();
    // draws the epistasis value for each PW postion pairs
    std::vector<double> drawEpistasis();
    // computes the index within the vector representation of a symmetric matrix
    unsigned int getMatrixVectorIndex(unsigned int i, unsigned int j);

public:

    const std::vector<double> &getKd();

    const double &getKd(const unsigned int i);

    const std::vector<double> &getEpistasis();

    const double &getEpistasis(const unsigned int i, const unsigned int j);

    static FunctionalSequence& get_instance();

    void writeKdsToFile(const std::string& filename);
    void writeEpistasisToFile(const std::string& filename);

    FunctionalSequence(FunctionalSequence const&) = delete;

    FunctionalSequence(FunctionalSequence&&) = delete;

};


#endif //DCABENCHMARK_FUNCTIONALSEQUENCE_HPP
