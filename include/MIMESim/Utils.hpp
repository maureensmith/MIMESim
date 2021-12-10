//
//  Utils.hpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 22.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#ifndef Utils_hpp
#define Utils_hpp

#include <filesystem>
#include <map>
#include <vector>

namespace fs = std::filesystem;

namespace utils
{
    typedef std::map<unsigned int, unsigned int> mutatedPositions;

    // TODO: refactoring -> besesr zu DCA und in DCABenchmark benutzen
    unsigned long nChoosek(const unsigned n, const unsigned k);
    std::vector<unsigned int> getBinaryRange(unsigned int maxRange, unsigned int L);
    std::vector<unsigned int> getMultinomialRange(unsigned int maxRange, unsigned int L, unsigned int q);
    // TODO QUESTION ist in Species vorhanden.... umschiften, ebenso wie mutPosToIndex?
    mutatedPositions specIdxToMutPos(const unsigned long specIdx, const unsigned int L, const unsigned int numSymbols,
                                     const std::vector<unsigned int> &nMutRange);
    /**
     * Read in parameters from a given parameter file in the given result directory. If there is no file, use default
     * parameters
     * TODO: entweder result ordner angeben als Muss, wo ggf die parameter liste drin ist
     */
    void readParameters(const fs::path &outputPath);

    /**
     * TODO noch abfragen, dass nur bei "" in den cout geschrieben werden soll?
     * Write parameters into a parameter file in the given output path.
     * If the path does not exist, the parameters are printed into cout
     */
    void writeParameters(const fs::path &outputPath);

    /**
     * Write parameters into cout
     */
    void writeParameters();
}

#endif /* Utils_hpp */
