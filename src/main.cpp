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

    fs::path outputPath("../results");
    if(argc > 1) {
        outputPath = argv[1];
    }

    if(!fs::exists(outputPath)){
        //fs::path outpath(outputPath);
        fs::create_directory(outputPath);
        std::cout << "Create output directory " << fs::canonical(outputPath) << std::endl;
    } else {
        std::cout << "Using output directory " << fs::canonical(outputPath) << std::endl;
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

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Create species *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
   // Create M species, the map contains the counts for each sampled sequence id
    species::species_map species_vec = species::drawSpeciesIds();

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";


//TODO: diese Art der Abfrage in die Tests packen
//    std::cout << "wt species count + freq. " << species_vec.at(1).getCount() << " " << species_vec.at(1).getFreq() << std::endl;
//    std::cout << "mut species count + freq. " << species_vec.at(2).getCount() << " " << species_vec.at(2).getFreq() << std::endl;
//    std::cout << "mut species count + freq. " << species_vec.at(3).getCount() << " " << species_vec.at(3).getFreq() << std::endl;
//    //std::cout << "mut bound unbound freq. " << species_vec.at(20877).getCount() << " " << species_vec.at(20877).getFreq() << std::endl;

    std::cout << "****** Create unmutated wild type library  *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    //auto wtSpecId_map = species::drawWildtypeErrors();
    //std::vector<species::Species> wtSpecies_vec;
    //The "control expereriment / wild type library" contains only wildtype sequences
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
    // create reference (here only As, Cs Gs and Us are mutations)
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
        //TODO workaround, anders lösen
        auto read = spec.second.createRead();
        counter_bound_1d.count(read, times_bound);
        counter_unbound_1d.count(read, times_unbound);

        counter_bound_2d.count(read, times_bound);
        counter_unbound_2d.count(read, times_unbound);

        //if (spec.first < 100) {
            //std::cout << "id " << spec.first << " countTot " << spec.second.getCount() << " bound " << times_bound
            //      << " unbound " << times_unbound << std::endl;
           // std::cout << " id: " << spec.first << " " << spec.second.getNumMut() << " " << spec.second.getCount() << " " <<  double(times_unbound)/spec.second.getCount()
            //          << " " << double(times_bound)/spec.second.getCount()  <<std::endl;
        //}

//        for(auto mut1:spec.second.getMutatedPositions()) {
//
//            //.... substract the ones where position 1 is mutated, and count as wt mut
//            for(unsigned mut2=1; mut2 <= cons.L; ++mut2) {
//                if (mut1.getPosition() < mut2) {
//                    counter_bound_2d.count(mut1, mut2, 'A', 'A', -times_bound);
//                    counter_bound_2d.count(mut1, mut2, 'C', 'A', times_bound);
//                    counter_unbound_2d.count(mut1, mut2, 'A', 'A', -times_unbound);
//                    counter_unbound_2d.count(mut1, mut2, 'C', 'A', times_unbound);
//                } else if(mut1 > mut2) {
//                    counter_bound_2d.count(mut2, mut1, 'A', 'A', -times_bound);
//                    counter_bound_2d.count(mut2, mut1, 'A', 'C', times_bound);
//                    counter_unbound_2d.count(mut2, mut1, 'A', 'A', -times_unbound);
//                    counter_unbound_2d.count(mut2, mut1, 'A', 'C', times_unbound);
//                }
//            }
//            //... substruct again from the wt mut, and at for mut mut... as this vector is very short anyway (max 4 or 5)
//            for(auto pos2:spec.second.getMutatedPositions()) {
//                if (mut1 < pos2) {
//                    counter_bound_2d.count(mut1, pos2, 'C', 'A', -times_bound);
//                    counter_bound_2d.count(mut1, pos2, 'A', 'C', -times_bound);
//                    counter_bound_2d.count(mut1, pos2, 'C', 'C', times_bound);
//                    counter_unbound_2d.count(mut1, pos2, 'C', 'A', -times_unbound);
//                    counter_unbound_2d.count(mut1, pos2, 'A', 'C', -times_unbound);
//                    counter_unbound_2d.count(mut1, pos2, 'C', 'C', times_unbound);
//                }
//
//                //TODO: hier wird doppelt gezählt, da ich 2 mal durch alle mutierten positionen gehe -> WEG
////                else if(mut1 > pos2) {
////                    counter_bound_2d.count(pos2, mut1, 'A', 'C', -times_bound);
////                    counter_bound_2d.count(pos2, mut1, 'C', 'C', times_bound);
////                    counter_unbound_2d.count(pos2, mut1, 'A', 'C', -times_unbound);
////                    counter_unbound_2d.count(pos2, mut1, 'C', 'C', times_unbound);
////                }
//            }
//
//        }

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

        //TODO workaround, ändern?
        auto read = spec.second.createRead();

        counter_bound_1d_wt.count(read , times_bound);
        counter_unbound_1d_wt.count(read, times_unbound);

        counter_bound_2d_wt.count(read, times_bound);
        counter_unbound_2d_wt.count(read, times_unbound);

//        if (spec.first < 100) {
//        std::cout << "id " << spec.first << " countTot " << spec.second.getCount() << " bound " << times_bound
//              << " unbound " << times_unbound << std::endl;
//            std::cout << " bound mut " << spec.second.getMutCountBound() <<  " bound err " << spec.second.getErrorCountBound() << std::endl;
//            std::cout << " unbound mut " << spec.second.getMutCountUnbound() <<  " unbound err " << spec.second.getErrorCountUnbound() << std::endl;
//        }

//        for(auto pos1:spec.second.getMutatedPositions()) {
//
//            for(unsigned pos2=1; pos2<=cons.L; ++pos2) {
//                if (pos1 < pos2) {
//                    counter_bound_2d_wt.count(pos1, pos2, 'A', 'A', -times_bound);
//                    counter_bound_2d_wt.count(pos1, pos2, 'C', 'A', times_bound);
//                    counter_unbound_2d_wt.count(pos1, pos2, 'A', 'A', -times_unbound);
//                    counter_unbound_2d_wt.count(pos1, pos2, 'C', 'A', times_unbound);
//                } else if(pos1 > pos2) {
//                    counter_bound_2d_wt.count(pos2, pos1, 'A', 'A', -times_bound);
//                    counter_bound_2d_wt.count(pos2, pos1, 'A', 'C', times_bound);
//                    counter_unbound_2d_wt.count(pos2, pos1, 'A', 'A', -times_unbound);
//                    counter_unbound_2d_wt.count(pos2, pos1, 'A', 'C', times_unbound);
//                }
//            }
//            for(auto pos2:spec.second.getMutatedPositions()) {
//                if (pos1 < pos2) {
//                    counter_bound_2d_wt.count(pos1, pos2, 'C', 'A', -times_bound);
//                    counter_bound_2d_wt.count(pos1, pos2, 'A', 'C', -times_bound);
//                    counter_bound_2d_wt.count(pos1, pos2, 'C', 'C', times_bound);
//                    counter_unbound_2d_wt.count(pos1, pos2, 'C', 'A', -times_unbound);
//                    counter_unbound_2d_wt.count(pos1, pos2, 'A', 'C', -times_unbound);
//                    counter_unbound_2d_wt.count(pos1, pos2, 'C', 'C', times_unbound);
//                }
//
//                //TODO: s.o.
////                else if(pos1 > pos2) {
////                    counter_bound_2d_wt.count(pos2, pos1, 'A', 'C', -times_bound);
////                    counter_bound_2d_wt.count(pos2, pos1, 'C', 'C', times_bound);
////                    counter_unbound_2d_wt.count(pos2, pos1, 'A', 'C', -times_unbound);
////                    counter_unbound_2d_wt.count(pos2, pos1, 'C', 'C', times_unbound);
////                }
//            }
//
//        }

    }

//     t2 = std::chrono::high_resolution_clock::now();
//    duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
//    std::cout <<"Time with L=" << L << " wtSpecVec " << wtSpecies_vec.size() << " " <<duration <<std::endl;

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";


    std::cout << "****** Write output to file *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    if(!fs::exists(outputPath));
        fs::create_directory(outputPath);

    //create subdirectories for the single and double mutant counts
    fs::create_directory(outputPath / "2d");
    fs::create_directory(outputPath / "1d");
    //create for each barcode a textfile and write the counts into the files
    counter_bound_1d.write_to_file(outputPath /"1d" / (mut_bound_id+".txt"));
    counter_unbound_1d.write_to_file(outputPath / "1d" / (mut_unbound_id+".txt"));
    counter_bound_2d.write_to_file(outputPath / "2d" / (mut_bound_id+".txt"));
    counter_unbound_2d.write_to_file(outputPath / "2d" / (mut_unbound_id+".txt"));

    counter_bound_1d_wt.write_to_file(outputPath  / "1d" / (wt_bound_id+".txt"));
    counter_unbound_1d_wt.write_to_file(outputPath / "1d" / (wt_unbound_id+".txt"));
    counter_bound_2d_wt.write_to_file(outputPath / "2d" / (wt_bound_id+".txt"));
    counter_unbound_2d_wt.write_to_file(outputPath / "2d" / (wt_unbound_id+".txt"));

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Write true values to files *******" << std::endl;
    effects.writeEpistasisToFile(outputPath / "pairwise_epistasis.txt");
    effects.writeKdsToFile(outputPath / "single_kds.txt");

    //write pairwise effects: species ids after 1(0 mut), and 2-L+1 (1 mut)
    std::ofstream outfile(outputPath / "pairwise_kds.txt");
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

    utils::writeParameters(cons.OUTPUT_DIR);

    end = std::chrono::high_resolution_clock::now();
    diff = end-start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Done *******" << std::endl;


    return 0;
}



