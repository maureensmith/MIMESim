//
//  main.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 03.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <chrono>

#include "Constants.hpp"
#include "Species.hpp"
#include "Utils.hpp"
#include "Mutation.hpp"
#include "FunctionalSequence.hpp"
#include "BindingCompetition.hpp"

#include "sam2counts/reference.hpp"
#include "sam2counts/nucleobase.hpp"
#include "sam2counts/count.hpp"

namespace fs = std::filesystem;

int main(int argc, const char * argv[]) {

    std::cout << "****** Set up constants *******" << std::endl;
    //measuring realtime duration (see std::clock for cpu time)
    auto start = std::chrono::high_resolution_clock::now();

    //TODO besser Lösung für .. finden
    std::string outputPath("../results");
    if(argc > 1) {
        outputPath = argv[1];
    }

    if(!fs::exists(outputPath)){
        fs::path outpath(outputPath);
        fs::create_directory(outpath);
        std::cout << "Create output directory " << absolute(outpath).string() << std::endl;
    } else {
        std::cout << "Using output directory " << outputPath << std::endl;
    }

    utils::readParameters(outputPath);

    // get the newly created instance of the constants
    auto& cons = constants::Constants::get_instance();

    // The 4 output files are saved with their ids where wild_type_bound = firstId, wild_type_unbound = firstId+1,
            // mut_bound = firstId+2, mut_unbound = firstId+3
    int firstId = 1;
    std::string wt_bound_id(std::to_string(firstId));
    std::string wt_unbound_id(std::to_string(firstId+1));
    std::string mut_bound_id(std::to_string(firstId+2));
    std::string mut_unbound_id(std::to_string(firstId+3));

//TODO in einen Test mit einbringen
//    std::cout << "Maximal number of mutations: " << cons.MAX_MUT << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Sample mutational effects *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    // Create Ground Truth: Effects of each mutated position and epistatic effects and sequencing noise
    FunctionalSequence& effects = FunctionalSequence::get_instance();
    std::vector<Mutation> mutationsPerPos_vec;
    mutationsPerPos_vec.reserve(cons.L);
    for(int i = 1; i<=cons.L; ++i) {
        //create instance for each
        mutationsPerPos_vec.emplace_back(i, effects.getKd(i));
    }

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Create species *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // Create M species
    auto specId_map = species::drawSpeciesIds();
    //std::vector<species::Species> species_vec;
    species::species_map species_vec;
    //species_vec.reserve(specId_map.size());
    for(auto it = specId_map.begin(); it != specId_map.end(); ++it) {
        // creates object directly in vector (instead of creating & moving), calling constr with given parameter
                //TODO: emplace_hint
        auto currentObj = species_vec.emplace(it->first, it->first);
        // first is a pointer to just constructed pair, second is the species object to call the methods
        currentObj.first->second.setCount(it->second);
        currentObj.first->second.computeSpeciesKd();
        //species_vec.emplace_back(it->first);
        //species_vec.back().setCount(it->second);
        ////species_vec.back().setErrors(it->second[1]);
        ////species_vec.back().setFreq(it->second[0]/double(cons.M));
        //species_vec.back().computeSpeciesKd();
        //std::cout << "id " << it->first << " nummut " << species_vec.at(it->first).getNumMut() << " counts " << it->second  << std::endl;
    }
    //TODO: diese Art der Abfrage in die Tests packen
//    std::cout << "wt species count + freq. " << species_vec.at(1).getCount() << " " << species_vec.at(1).getFreq() << std::endl;
//    std::cout << "mut species count + freq. " << species_vec.at(2).getCount() << " " << species_vec.at(2).getFreq() << std::endl;
//    std::cout << "mut species count + freq. " << species_vec.at(3).getCount() << " " << species_vec.at(3).getFreq() << std::endl;
//    //std::cout << "mut bound unbound freq. " << species_vec.at(20877).getCount() << " " << species_vec.at(20877).getFreq() << std::endl;

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Create unmutated wild type library  *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    //auto wtSpecId_map = species::drawWildtypeErrors();
    //std::vector<species::Species> wtSpecies_vec;
    species::species_map wtSpecies_vec;
    auto currentObj = wtSpecies_vec.emplace(1, 1);
    currentObj.first->second.setCount(cons.M);
    currentObj.first->second.computeSpeciesKd();

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";


    std::cout << "****** Solve ODE to infer bound and unbound fraction of species *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    //TODO: Umbau nach counts
    //std::valarray<double> f_bound_tot;
    //std::valarray<double> f_unbound_tot;
    std::valarray<unsigned int> S_bound(species_vec.size());
    std::valarray<unsigned int> S_unbound(species_vec.size());

    //std::valarray<int> f_bound_tot(species_vec.size());
    //std::valarray<int> f_unbound_tot(species_vec.size());
    // set up the ODE (binding competition) and solve it to get the bound and unbound fractions (from the total amount M) in equilibrium
    UnboundProtein f(species_vec);
    f.solve(S_bound, S_unbound);

    //std::cout << "wt bound unbound freq. " << f_bound_tot[0] << " " << f_unbound_tot[0] << std::endl;
    //std::cout << "mut bound unbound freq. " << f_bound_tot[1] << " " << f_unbound_tot[1] << std::endl;
    //std::cout << "mut bound unbound freq. " << f_bound_tot[2] << " " << f_unbound_tot[2] << std::endl;


    //*************solve for wt
    //TODO: Umbau nach counts
    //std::valarray<double> f_bound_tot_wt;
    //std::valarray<double> f_unbound_tot_wt;
    std::valarray<unsigned int> S_bound_wt(wtSpecies_vec.size());
    std::valarray<unsigned int> S_unbound_wt(wtSpecies_vec.size());
    // set up the ODE (binding competition) and solve it to get the bound and unbound fractions (from the total amount M) in equilibrium
    UnboundProtein f_wt(wtSpecies_vec);
    f_wt.solve(S_bound_wt, S_unbound_wt);

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Add noise *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    //Carefull: The map is extended by species that occur only because of sequencing error, hence the length of S_bound and
    //S_unbound dont fit any more with the length of the map
    //TODO: Umbau nach counts
    species::addCountsWithError(S_bound, S_unbound, species_vec);

    species::addCountsWithError(S_bound_wt, S_unbound_wt, wtSpecies_vec);

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Count Mut Library *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // create reference (here only As, Cs are mutations)
    ref::reference ref;
    ref::ref_map wt_read;
    wt_read.reserve(cons.L);
    for(int i=1; i <=cons.L; ++i) {
        ref.add(nucleotid::nucleobase{1});
        wt_read.add({i, nucleotid::nucleobase{'A'}});
    }
    //count for mutation library
    count::counter_1 counter_bound_1d{ref};
    count::counter_1 counter_unbound_1d{ref};

    //Because the majority is wildtype, count all als wildtype....
    count::counter_2 counter_bound_2d{ref};
    //TODO: Umbau nach counts
    //counter_bound_2d.count(wt_read, round(cons.M*f_bound_tot.sum()));
    counter_bound_2d.count(wt_read, S_bound.sum());
    count::counter_2 counter_unbound_2d{ref};
    //TODO: Umbau nach counts
    //counter_unbound_2d.count(wt_read, round(cons.M*f_unbound_tot.sum()));
    counter_unbound_2d.count(wt_read, S_unbound.sum());
    //count::counter_3 counter_bound{ref};
    //count::counter_3 counter_unbound{ref};


    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();



    //for each species: set up read, and let it count the nucleotides
     int i = 0;
    for(auto& spec : species_vec) {
//        std::cout << "bound "<< spec.second.getMutCountBound()<< std::endl;
//        std::cout << "bound error " << spec.second.getErrorCountBound() << std::endl;
//        std::cout << "unbound "<< spec.second.getMutCountUnbound() << std::endl;
//        std::cout << "unbound error "<< spec.second.getErrorCountUnbound() << std::endl;
        //TODO auch noch schneller machen

        const int times_bound = spec.second.getMutCountBound() + spec.second.getErrorCountBound();
        const int times_unbound = spec.second.getMutCountUnbound() + spec.second.getErrorCountUnbound();
        counter_bound_1d.count(spec.second.getRead(), times_bound);
        counter_unbound_1d.count(spec.second.getRead(), times_unbound);

        //if (spec.first < 100) {
            //std::cout << "id " << spec.first << " countTot " << spec.second.getCount() << " bound " << times_bound
            //      << " unbound " << times_unbound << std::endl;
           // std::cout << " id: " << spec.first << " " << spec.second.getNumMut() << " " << spec.second.getCount() << " " <<  double(times_unbound)/spec.second.getCount()
            //          << " " << double(times_bound)/spec.second.getCount()  <<std::endl;
        //}

        for(auto pos1:spec.second.getMutatedPositions()) {

            //.... substract the ones where position 1 is mutated, and count as wt mut
            for(unsigned pos2=1; pos2<=cons.L; ++pos2) {
                if (pos1 < pos2) {
                    counter_bound_2d.count(pos1, pos2, 'A', 'A', -times_bound);
                    counter_bound_2d.count(pos1, pos2, 'C', 'A', times_bound);
                    counter_unbound_2d.count(pos1, pos2, 'A', 'A', -times_unbound);
                    counter_unbound_2d.count(pos1, pos2, 'C', 'A', times_unbound);
                } else if(pos1 > pos2) {
                    counter_bound_2d.count(pos2, pos1, 'A', 'A', -times_bound);
                    counter_bound_2d.count(pos2, pos1, 'A', 'C', times_bound);
                    counter_unbound_2d.count(pos2, pos1, 'A', 'A', -times_unbound);
                    counter_unbound_2d.count(pos2, pos1, 'A', 'C', times_unbound);
                }
            }
            //... substruct again from the wt mut, and at for mut mut... as this vector is very short anyway (max 4 or 5)
            for(auto pos2:spec.second.getMutatedPositions()) {
                if (pos1 < pos2) {
                    counter_bound_2d.count(pos1, pos2, 'C', 'A', -times_bound);
                    counter_bound_2d.count(pos1, pos2, 'A', 'C', -times_bound);
                    counter_bound_2d.count(pos1, pos2, 'C', 'C', times_bound);
                    counter_unbound_2d.count(pos1, pos2, 'C', 'A', -times_unbound);
                    counter_unbound_2d.count(pos1, pos2, 'A', 'C', -times_unbound);
                    counter_unbound_2d.count(pos1, pos2, 'C', 'C', times_unbound);
                }

                //TODO: hier wird doppelt gezählt, da ich 2 mal durch alle mutierten positionen gehe -> WEG
//                else if(pos1 > pos2) {
//                    counter_bound_2d.count(pos2, pos1, 'A', 'C', -times_bound);
//                    counter_bound_2d.count(pos2, pos1, 'C', 'C', times_bound);
//                    counter_unbound_2d.count(pos2, pos1, 'A', 'C', -times_unbound);
//                    counter_unbound_2d.count(pos2, pos1, 'C', 'C', times_unbound);
//                }
            }

        }

        //counter_bound_2d.count(spec.second.getRead(), spec.second.getMutCountBound() + spec.second.getErrorCountBound());
        //counter_unbound_2d.count(spec.second.getRead(), spec.second.getMutCountUnbound() + spec.second.getErrorCountUnbound());
        //++i;
        //if(i == 250000)
            //break;
    }

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

//    std::cout << std::endl;
//    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();

   // std::cout <<"Time with L=" << L << " specVec " << species_vec.size() << " " <<duration <<std::endl;

    std::cout << "****** Count Wt Library *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    //count for wt
    count::counter_1 counter_bound_1d_wt{ref};
    count::counter_1 counter_unbound_1d_wt{ref};

    count::counter_2 counter_bound_2d_wt{ref};
    //TODO: Umbau nach counts
    //counter_bound_2d_wt.count(wt_read, round(cons.M*f_bound_tot_wt.sum()));
    counter_bound_2d_wt.count(wt_read, S_bound_wt.sum());
    count::counter_2 counter_unbound_2d_wt{ref};
    //TODO: Umbau nach counts
    //counter_unbound_2d_wt.count(wt_read, round(cons.M*f_unbound_tot_wt.sum()));
    counter_unbound_2d_wt.count(wt_read, S_unbound_wt.sum());

    t1 = std::chrono::high_resolution_clock::now();

    for(auto& spec : wtSpecies_vec) {
        //TODO auch noch schneller machen

        const int times_bound = spec.second.getMutCountBound() + spec.second.getErrorCountBound();
        const int times_unbound = spec.second.getMutCountUnbound() + spec.second.getErrorCountUnbound();

        counter_bound_1d_wt.count(spec.second.getRead(), times_bound);
        counter_unbound_1d_wt.count(spec.second.getRead(), times_unbound);

//        if (spec.first < 100) {
//        std::cout << "id " << spec.first << " countTot " << spec.second.getCount() << " bound " << times_bound
//              << " unbound " << times_unbound << std::endl;
//            std::cout << " bound mut " << spec.second.getMutCountBound() <<  " bound err " << spec.second.getErrorCountBound() << std::endl;
//            std::cout << " unbound mut " << spec.second.getMutCountUnbound() <<  " unbound err " << spec.second.getErrorCountUnbound() << std::endl;
//        }

        for(auto pos1:spec.second.getMutatedPositions()) {

            for(unsigned pos2=1; pos2<=cons.L; ++pos2) {
                if (pos1 < pos2) {
                    counter_bound_2d_wt.count(pos1, pos2, 'A', 'A', -times_bound);
                    counter_bound_2d_wt.count(pos1, pos2, 'C', 'A', times_bound);
                    counter_unbound_2d_wt.count(pos1, pos2, 'A', 'A', -times_unbound);
                    counter_unbound_2d_wt.count(pos1, pos2, 'C', 'A', times_unbound);
                } else if(pos1 > pos2) {
                    counter_bound_2d_wt.count(pos2, pos1, 'A', 'A', -times_bound);
                    counter_bound_2d_wt.count(pos2, pos1, 'A', 'C', times_bound);
                    counter_unbound_2d_wt.count(pos2, pos1, 'A', 'A', -times_unbound);
                    counter_unbound_2d_wt.count(pos2, pos1, 'A', 'C', times_unbound);
                }
            }
            for(auto pos2:spec.second.getMutatedPositions()) {
                if (pos1 < pos2) {
                    counter_bound_2d_wt.count(pos1, pos2, 'C', 'A', -times_bound);
                    counter_bound_2d_wt.count(pos1, pos2, 'A', 'C', -times_bound);
                    counter_bound_2d_wt.count(pos1, pos2, 'C', 'C', times_bound);
                    counter_unbound_2d_wt.count(pos1, pos2, 'C', 'A', -times_unbound);
                    counter_unbound_2d_wt.count(pos1, pos2, 'A', 'C', -times_unbound);
                    counter_unbound_2d_wt.count(pos1, pos2, 'C', 'C', times_unbound);
                }

                //TODO: s.o.
//                else if(pos1 > pos2) {
//                    counter_bound_2d_wt.count(pos2, pos1, 'A', 'C', -times_bound);
//                    counter_bound_2d_wt.count(pos2, pos1, 'C', 'C', times_bound);
//                    counter_unbound_2d_wt.count(pos2, pos1, 'A', 'C', -times_unbound);
//                    counter_unbound_2d_wt.count(pos2, pos1, 'C', 'C', times_unbound);
//                }
            }

        }

    }

//     t2 = std::chrono::high_resolution_clock::now();
//    duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
//    std::cout <<"Time with L=" << L << " wtSpecVec " << wtSpecies_vec.size() << " " <<duration <<std::endl;

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";


    std::cout << "****** Write output to file *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    std::filesystem::path outPath =outputPath;
    if(!fs::exists(outPath));
        fs::create_directory(outPath);
    fs::create_directory(outputPath+"/2d");
    fs::create_directory(outputPath+ "/1d");
    counter_bound_1d.write_to_file(outputPath +"/1d/"+mut_bound_id+".txt");
    counter_unbound_1d.write_to_file(outputPath +"/1d/"+mut_unbound_id+".txt");
    counter_bound_2d.write_to_file(outputPath +"/2d/"+mut_bound_id+".txt");
    counter_unbound_2d.write_to_file(outputPath +"/2d/"+mut_unbound_id+".txt");

    counter_bound_1d_wt.write_to_file(outputPath +"/1d/"+wt_bound_id+".txt");
    counter_unbound_1d_wt.write_to_file(outputPath +"/1d/"+wt_unbound_id+".txt");
    counter_bound_2d_wt.write_to_file(outputPath +"/2d/"+wt_bound_id+".txt");
    counter_unbound_2d_wt.write_to_file(outputPath +"/2d/"+wt_unbound_id+".txt");

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Write true values to files *******" << std::endl;
    effects.writeEpistasisToFile(outputPath+"/pairwise_epistasis.txt");
    effects.writeKdsToFile(outputPath+"/single_kds.txt");

    //write pairwise effects: species ids after 1(0 mut), and 2-L+1 (1 mut)
    std::ofstream outfile(outputPath+"/pairwise_kds.txt");
    if (outfile.good())
    {

        for(int i=cons.NMUT_RANGE[1]+1; i<=cons.NMUT_RANGE[2];++i) {
            outfile << species_vec.at(i).getKd() << '\n';
        }
        outfile.close();
    }

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Write parameter into File *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    utils::writeParameters();

    std::ofstream paraOutfile(outputPath+"/parameter.txt");

    if (outfile.good())
    {
        paraOutfile << "### paratemers regardin kd sampling ###\n";
        paraOutfile << "kd_wt\t" <<  cons.KD_WT << '\n';
        paraOutfile << "p_effect\t" <<  cons.P_EFFECT << '\n';
        paraOutfile << "p_epistasis\t" <<  cons.P_EPISTASIS << '\n';
        paraOutfile << "### paramters regardin kd sampling ###\n";
        paraOutfile << "L\t" <<  cons.L << '\n';
        paraOutfile << "q\t" <<  cons.Q << '\n';
        paraOutfile << "M\t" <<  cons.M << '\n';
        paraOutfile << "p_mut\t" <<  cons.P_MUT << '\n';
        paraOutfile << "p_error\t" <<  cons.P_ERR << '\n';
        paraOutfile.close();
    }

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Done *******" << std::endl;


    return 0;
}



