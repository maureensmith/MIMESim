# MIMESim

1. Simulate sequence variants in pool of a wild type sequence (random mutagenesis)
2. Randomly assign fitness values for each position-wise mutation and pairwise epistatic effects
3. Derive species frequencies in equilibrium
4. Add statistical noise
5. Output positionwise mutant counts and the groundtruth

*Detailed explanation will be added*

## Install

Go into the MIMESim directory and create a build directory and cmake and make from there.

```
cd /path/to/MIMESim
mkdir build
cd build
cmake ..
make
```
This will create the binary file *.../MIMESim/build/bin/MIMESim_prog*

## Run

Running the program with an optional parameter, giving the path to the result directory. 
If no parameter is given, the output will be written to the newly created directory *../results*.

```
/path/to/MIMESim_prog [result_directory]
```

If no parameter file is present in the given result directory (or no directory was given at all), the programm will run the simulation with default paraters *(tba)*.
After the run it will create the parameter file in the result directory together with the other output. 

If you want to set the simulation parameters, create the *parameters.txt* file in the desired result directory and add the following parameter names and values (tab separated):

| parameter       | type          | description  |
| :-------------  |:-------------:| :------------|
| kd_wt           | (float)       |   (absolute Kd for the wildtype sequence (default: 1))|
| p_effect        | (float)       |   probability per position to have an effect on function when mutated (default: 0.5) |
| p_epistasis     | (float)       |   probability for each pair of mutations to have epistatic effect (default: 0.3) |
| L               | (int)         |   Length of the sequence (default: 50) |
| q               | (int)         |   (number of possible symbols per position (default: 2))  |
| M               | (int)         |   (size of initial population (default: 12000000)  |
| p_mut           | (float)       |   probability for each position of being mutated (default: 0.1) |
| p_error         | (float)       |   probability of introducing an error (default: 0.01) |

**Comments:**
1. Right now, kd_wt is always 1 and also M is always 12 mio, but this can be easily adjusted in the code.
2. It is only working for q=2 atm.
3. B_tot = total protein concentration is not in the parameter setting right now.
The amount of protein is given in relation to the amount of sequences, and is set to 2.0


