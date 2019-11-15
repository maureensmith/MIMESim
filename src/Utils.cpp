//
//  Utils.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 22.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#include "Utils.hpp"
#include "Constants.hpp"
#include <iostream>
#include <math.h>
#include <numeric>
#include <filesystem>
#include <iostream>
#include <fstream>

namespace fs = std::filesystem;

//TODO: Utils in DCA packen, da es nicht Benchmark spezifisch ist

namespace utils {

    unsigned long nChoosek(unsigned n, unsigned k)
    {
        if (k > n) return 0;
        if (k * 2 > n) k = n-k;
        if (k == 0) return 1;

        unsigned long result = n;
        for( int i = 2; i <= k; ++i ) {
            result *= (n-i+1);
            result /= i;
        }
        return result;
    }

    std::vector<unsigned int> getBinaryRange(unsigned int maxRange, unsigned int L) {
        //compute the number of possible sequence for 0..MAX_MUT mutations, the cumulative sum gives the id range for each number of mutations
        std::vector<unsigned int> nMutRange(maxRange+1);
        nMutRange[0] = 1;
        for(unsigned int i = 1; i<=maxRange; ++i) {
            nMutRange[i] = nMutRange[i - 1] + utils::nChoosek(L, i);
        }
        return nMutRange;
    }

    std::vector<unsigned int> getMultinomialRange(unsigned int maxRange, unsigned int L, unsigned int q) {
        //compute the number of possible sequence for 0..MAX_MUT mutations, the cumulative sum gives the id range for each number of mutations
        std::vector<unsigned int> nMutRange(maxRange+1);
        nMutRange[0] = 1;
        for(unsigned int i = 1; i<=maxRange; ++i) {
            nMutRange[i] = nMutRange[i - 1] + (utils::nChoosek(L, i) * std::pow(q,i));
        }
        return nMutRange;
    }


    mutatedPositions specIdxToMutPos(const unsigned long specIdx, const unsigned int L, const unsigned int numSymbols,
                                     const std::vector<unsigned int> &nMutRange) {
        mutatedPositions mutPos{};
        // check if Id is valid
        if(specIdx <= nMutRange.back() && specIdx >0) {
            // number of mutations in sequence
            unsigned int numMut = std::lower_bound(nMutRange.begin(), nMutRange.end(), specIdx) - nMutRange.begin();
            //std::cout << "id " << specIdx << " numMut " << numMut << std::endl;
            //collect the mutated positions and the mutation
            //std::map<unsigned int, unsigned int> pos;
            //std::vector<unsigned int> pos(numMut);
            //std::vector<unsigned int> mut(numMut);
            //id 1 = no mutations: don't start to compute anything
            if(numMut != 0) {
                unsigned int Lact = L;
                unsigned int numMutAct = numMut;
                //get the id within the range of number of mutations (substract the ids for the sequences with less mutations)
                unsigned long idAct = specIdx - nMutRange[numMut-1];
                //determine the each mutations position seen from the mutations position before...
                unsigned int m = 0;
                do{
                    //for each possible positions within the length the actual mutation covers a range of ids depending on the residual mutations to follow

                    //cumSumRange <- cumsum(choose(Lact-seq(Lact-(numMutAct-1)), numMutAct-1))
                    std::vector<long> cumSumRange(Lact-(numMutAct-1));
                    auto n = Lact - 1;
                    for(unsigned int i =0; i < cumSumRange.size(); ++i) {
                        //the id range for each position when it is mutated (if pos 1 is mutated L-1 positions remain to
                        // contain the mut-1 remaining mutations, each time there are q^nmut possibilies (it does not change in every step)
                        cumSumRange[i] = utils::nChoosek(n, numMutAct - 1) * std::pow(numSymbols,numMut);
                        --n;
                        //std::cout << "cumSumRange von " << i << " " << cumSumRange[i] << std::endl;
                    }

                    std::partial_sum(cumSumRange.begin(), cumSumRange.end(), cumSumRange.begin());

                    // find the id within the ranges and get the index (=position seen from the last mutated position)
                    int pos = std::lower_bound(cumSumRange.begin(), cumSumRange.end(), idAct) - cumSumRange.begin()+1;
                    // if one mutation, the symbol can be simply find by modulo q, with 2 mutations it has to be divided by q first
                            // (since one symbol of the first mutation can be present with q symbols of the second one, with 3 mutation divided by q^2 etc
                     int mut = (int((idAct - 1) / std::pow(numSymbols,numMutAct-1))  % numSymbols ) + 1;
                    //put the position - symbol pair into the map. add the last cummulative position (=position seen from the beginning of the sequence)
                    unsigned int prePos = 0;
                    if(mutPos.begin() != mutPos.end())
                        prePos = mutPos.rbegin()->first;
                    mutPos.insert({pos+ prePos, mut});
                    //std::cout << "mut " << m << ": " << pos <<" " << mut << std::endl;
                    // the residual length after the actual mutation
                    Lact = Lact - pos;
                    // the redsiudal number of mutations after the actual mutation
                    --numMutAct;
                    // the id within the residual length (-1 because the inidices start at 0, and another -1 because we substract the ids of the preceding mutation range
                    //changees for q symbols: the substracted ids have to be multiplied by q (since each position can have 3 values)
                    idAct = idAct - (pos==1? 0 : cumSumRange[pos-2]);
                    ++m;
                }while(m<numMut);
//                // .... and get the correct positions within the sequence with cumsum
//                std::partial_sum(pos.begin(), pos.end(), pos.begin());
//            for(auto p : pos)
//                std::cout << " " << p;
                //std::cout << std::endl;
            }

        } else {
            std::cerr << "ERROR: id is too big for length of sequence and number of mutations" << std::endl;
        }
        return mutPos;
    }

    // use only locally
    static void printParameters() {
        auto& cons = constants::Constants::get_instance();

        std::cout << std::endl;
        std::cout << "### paratemers regarding kd sampling ###\n";
        std::cout << "kd_wt\t" <<  cons.KD_WT << '\n';
        std::cout << "p_effect\t" <<  cons.P_EFFECT << '\n';
        std::cout << "p_epistasis\t" <<  cons.P_EPISTASIS << '\n';
        std::cout << "### paramters regarding sequence sampling ###\n";
        std::cout << "L\t" <<  cons.L << '\n';
        std::cout << "q\t" <<  cons.Q << '\n';
        std::cout << "M\t" <<  cons.M << '\n';
        std::cout << "p_mut\t" <<  cons.P_MUT << '\n';
        std::cout << "p_error\t" <<  cons.P_ERR << '\n';
        std::cout << std::endl;
    }


    void readParameters(const std::string &outputPath){

        // if output directory does not exist yet, create it
        if(!fs::exists(outputPath)){
            fs::path outpath(outputPath);
            fs::create_directory(outpath);
            std::cout << "Create output directory " << absolute(outpath).string() << std::endl;
        }


        // dafault parameters, in case no file is given, or paramater are not set in the file
        //TODO If L is big, the id range gets really high -> long statt int. if any error, look for to high ints
        //sequence length
        unsigned int L = 50;
        // symbols per position
        unsigned int q = 2;
        //mutation probability
        double p_mut = 0.1;
        double p_error = p_mut / 10.0;
        double p_effect = 0.5;
        double p_epistasis = 0.3;
        //unsigned int M = 12 * pow(10, 6);

        const fs::path paraFile = outputPath+"/parameter.txt";
        // if the output directory contains the parameter file, read it
        if(fs::exists(paraFile) && fs::is_regular_file(paraFile)) {
            std::ifstream infile(paraFile.string());

            if (infile.good()) {
                std::cout << "Reading file " << paraFile.string();
                std::string line;
                std::string param;
                std::string val;


                try {

                    while (std::getline(infile, line)) {

                        //stream through each line to read the parameter name and its value, seperated by tabular
                        std::istringstream lineSS(line);

                        std::getline(lineSS, param, '\t');
                        std::getline(lineSS, val, '\t');
                        //TODO ist erstmal fest auf 1.0 gesetzt und M auf 12 Mio
                        //if(param == "kd_wt")
                        //    kd_wt = std::stoi(val);
                        //if(param == "M")
                        //    M = std:stoi(val);
                        if (param == "L")
                            L = std::stoi(val);
                        if (param == "q")
                            q = std::stoi(val);
                        if (param == "p_mut")
                            p_mut = std::stod(val);
                        if (param == "p_error")
                            p_error = std::stod(val);
                        if (param == "p_effect")
                            p_effect = std::stod(val);
                        if (param == "p_epistasis")
                            p_epistasis = std::stod(val);
                    }
                    infile.close();
                    std::cout << " ... successful." << std::endl;
                }  catch (const std::invalid_argument& ia) {
                    //TODO anders mit umgehen?
                    std::cerr << "Error in parameter file. Invalid argument for parameter " << param << " ("<< ia.what() << ")"<< '\n';
                    std::cout << "Using default parameters." << std::endl;
                }
            }
        } else {
            //...if no config file is given, create the constants with the default parameter values
            std::cout << "No parameter file given. Using default parameters." << std::endl;
            //TODO so oder so noch die parameter ausgeben.
            //TODO test case für parameter einlesen
        }

        // Create constants which are used through out this test set
        constants::Constants& cons = constants::Constants::create_instance(L, q, p_mut, p_error, p_effect, p_epistasis);
        printParameters();
    }


    void writeParameters() {

        auto& cons = constants::Constants::get_instance();

        std::ofstream paraOutfile(cons.PARAMETER_FILE);

        if (paraOutfile.good())
        {
            paraOutfile << "### paratemers regarding kd sampling ###\n";
            paraOutfile << "kd_wt\t" <<  cons.KD_WT << '\n';
            paraOutfile << "p_effect\t" <<  cons.P_EFFECT << '\n';
            paraOutfile << "p_epistasis\t" <<  cons.P_EPISTASIS << '\n';
            paraOutfile << "### paramters regarding sequence sampling ###\n";
            paraOutfile << "L\t" <<  cons.L << '\n';
            paraOutfile << "q\t" <<  cons.Q << '\n';
            paraOutfile << "M\t" <<  cons.M << '\n';
            paraOutfile << "p_mut\t" <<  cons.P_MUT << '\n';
            paraOutfile << "p_error\t" <<  cons.P_ERR << '\n';
            paraOutfile.close();
        }
    }

}


