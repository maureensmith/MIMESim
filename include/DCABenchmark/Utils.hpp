//
//  Utils.hpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 22.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#ifndef Utils_hpp
#define Utils_hpp

#include <vector>
#include <map>

namespace utils {
    typedef std::map<unsigned int, unsigned int> mutatedPositions;

    //TODO: refactoring -> besesr zu DCA und in DCABenchmark benutzen
    unsigned long nChoosek( const unsigned n, const unsigned k);
    std::vector<unsigned int> getBinaryRange(unsigned int maxRange, unsigned int L);
    std::vector<unsigned int> getMultinomialRange(unsigned int maxRange, unsigned int L, unsigned int q);
    mutatedPositions specIdxToMutPos(const unsigned long specIdx, const unsigned int L, const unsigned int numSymbols,
                                     const std::vector<unsigned int> &nMutRange);
    /**
     * Read in parameters from a given parameter file in the given result directory. If there is no file, use default parameters
     * TODO: entweder result ordner angeben als Muss, wo ggf die parameter liste drin ist
     */
    void readParameters(const std::string &outputPath);

    /**
     * Write parameters into parameters.txt file in result directory
     */
    void writeParameters();
}

#endif /* Utils_hpp */
