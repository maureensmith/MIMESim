//
// Created by Smith, Maureen on 31.05.18.
//

#ifndef Mutation_hpp
#define Mutation_hpp


class Mutation {
private:
    const unsigned int position;
    const double kd;
    unsigned int count;
    double freq;
    double noise;

public:
    Mutation(const unsigned int p, const double kd);

    const unsigned int getPosition() const;

    const double getKd() const;

    unsigned int getCount() const;

    double getFreq() const;

    void setCount(unsigned int count);

    void setFreq(double freq);

    double getNoise() const;

    void setNoise(double noise);

};


#endif //DCABENCHMARK_MUTATION_H
