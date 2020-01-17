//
// Created by Smith, Maureen on 31.05.18.
//

#ifndef Mutation_hpp
#define Mutation_hpp


class Mutation {
private:
    unsigned int position; //const! not declared here, to be able to put it into a set
    //TODO auch als char? (not necessarily, weil ja auch beliebiges alphabet gelten kann.
    unsigned int symbol; //const! not declared here, to be able to put it into a set
    //TODO kann weg?
    //const double kd;
//    unsigned int count;
//    double freq;
//    //TODO noise count? statt double?
//    double noise;

public:
    Mutation(const unsigned int p, const unsigned int s);

    //position from 1 .. L
    const unsigned int getPosition() const;

    //mutation symbol from 0 .. Q-1 -1
    //TODO ändern? von 1 .. Q-1?
    const unsigned int getSymbol() const;

    //TODO weg?
    //const double getKd() const;

//    unsigned int getCount() const;
//
//    double getFreq() const;
//
//    void setCount(unsigned int count);
//
//    void setFreq(double freq);
//
//    double getNoise() const;
//
//    void setNoise(double noise);

    //to be able to build a set of Mutations (where the position is unique), the following operators have to be overloaded.
    bool operator< (const Mutation& mut) const;
};


#endif //DCABENCHMARK_MUTATION_H
