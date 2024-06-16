#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "particle.h"
#include "random.h"

using namespace std;
using namespace arma;

class Metropolis {


private:
  const int _ndim = 3;  // Dimensionality of the system
  bool _restart;        // Flag indicating if the simulation is restarted
  int _nblocks;         // Number of blocks for block averaging
  int _n_eq_blocks;     // Number of blocks for equilibration
  int _nsteps;          // Number of simulation steps in each block
  int _n_eq_steps;      // Number of steps for each equilibration block 
  int _nattempts;       // Number of attempted moves
  int _naccepted;       // Number of accepted moves
  int _n_iter;           // number of iteration metropolis
  double _temp, _beta;  // Temperature and inverse temperature
  double _temp_restarted;// Temperature after equilibration
  double _metropolis_step_0;  // max lenght of metropolis step 0
  double _metropolis_step_1;  // max lenght of metropolis step 1
  Random _rnd;          // Random number generator
  vec _fx, _fy, _fz;    // Forces on particles along x, y, and z directions
  double _trajx, _trajy, _trajz; //trajectories
  double _block_av;         // Block averages of properties
  double _global_av;        // Global averages of properties
  double _global_av2;       // Squared global averages of properties
  double _average;          // Average values of properties
  vec _measurement;      // Measured values of properties
  double _block_av_eq;         // Block averages of properties eq
  double _global_av_eq;        // Global averages of properties eq
  double _global_av2_eq;       // Squared global averages of properties eq
  double _average_eq;          // Average values of properties eq
  double _x0;
  double _y0;
  double _z0;
  double _x;
  double _y;
  double _z; 
  double _x0_new;
  double _y0_new;
  double _z0_new;    // starting point
  bool _measure_r;
  bool _print_r;
  bool _print_r_first;
  double _sigma_0;  //sigma gaussian for state 100
  double _sigma_1; //sigma gaussian for state 210

public: // Function declarations

  int get_nbl();                                       // Get the number of blocks
  int get_nsteps();                                   // Get the number of steps in each block
  int get_nbl_eq();                                   // Get the number of blocks eq
  int get_nsteps_eq();                              // Get the number of steps in each block eq
  void Initialize(int state);                                  // initialize system giving initial position and numbero of steps/blocks
  void initialize_properties(int state);
  void block_reset(int blk);  // Reset block averages
  void newMetropolis_algorithm(int state);
  void Metropolis_algorithm(double x0, double y0, double z0, int state);      // Start the algorithm
  void Metropolis_algorithm_Gauss(int state); // Start the algorithm, with normal transition probability
  int get_nattempts();                            // Get number of attemps
  int get_naccepted();                             // Get number of step accepted
  double mod_psi_0(double x, double y, double z);   // modulo quadro psi_100
  double mod_psi_1(double x, double y, double z);    // modulo quadro psi_210
  void r_averages(int blk, int state);                                // Calculate r_mean with block rule
  double get_initial_x();   
  double get_initial_y();
  double get_initial_z();                 // restituisce x0, y0 e z0
  double error(double acc, double acc2, int blk);     // errore blocking
  void set_nattempts();       // azzera n_attemps
  void set_naccepted();       // azzera n_accepted
  void set_r0();              // setta il punto di partenza a new x0, y0 e z0
  void print_traj(int state);
  void set_xyz(double x , double y, double z);
};



