#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "random.h"

using namespace std;
using namespace arma;

class Var_MC {
public:
    

    // Functions for computing various physical quantities
    double Gaussiana(double x, double mu, double sigma);
    double Kinetic(double x, double sigma, double mu);
    double Potenzial(double x);
    double Hamiltonian(double x, double sigma, double mu);
    double Psi(double x, double sigma, double mu);
    double mod_psi(double x, double sigma, double mu);

    // Metropolis algorithm
    void Metropolis_algorithm(double mu, double sigma);

    // Simulated Annealing
    void SimulatedAnnealing();

    // Initialization functions
    void initialize_properties();
    void Initialize();

    // Helper functions
    void print_traj();
    void block_reset(int blk);

    void set_temp(double temp);               // Set the temperature
    int get_nbl();
    int get_nsteps();
    double get_temp();
    int get_naccepted();
    int get_nattempts();
    int get_n_temps();
    double get_temp_up();
    double get_temp_low();
    double get_mu();
    double get_sigma();
    double get_x();

    void set_mu(double mu);
    void set_sigma(double sigma);
    void set_best_parameters(string filename);

    void measure();
    void averages(int blk);

    double error(double acc, double acc2, int blk);

    

    

private:
    // Parameters
    double _x;
    double _mu;
    double _sigma;
    double _metropolis_step;
    double _delta_opt;
    double _temp;
    double _t_low;
    double _t_up;
    double _n_temps;
    int _nblocks;
    int _nsteps;
    int _n_iter;
    int _nattempts;
    int _naccepted;
    int _SA_attempted;
    int _SA_accepted;
    bool _measure_penergy;
    bool _measure_kenergy;
    bool _measure_tenergy;
    int _temp_rep;
    int _nprop=0;

    

    // Random number generator
    Random _rnd;

    vec _block_av;         // Block averages of properties
    vec _global_av;        // Global averages of properties
    vec _global_av2;       // Squared global averages of properties
    vec _average;          // Average values of properties
    vec _measurement;      // Measured values of properties
    int _index_penergy, _index_kenergy, _index_tenergy;       // Indices for accessing energy-related properties in vec _measurement
};


