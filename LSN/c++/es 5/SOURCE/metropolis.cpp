#include <cmath>
#include <cstdlib>
#include <string>
#include "metropolis.h"
#include <armadillo>

using namespace std;
using namespace arma;



void Metropolis :: Initialize(int state){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);

  ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
  couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;
  couta.close();

  

  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat");
  string property;
  while ( !input.eof() ){
    input >> property;
    if (property=="X0"){
        input >> _x0;
        coutf << "X0" << _x0 << endl;
        cerr << "x0 ok" << endl;
        //_trajx.push_back(_x0);
    }
    else if (property=="Y0"){
        input >> _y0;
        coutf << "Y0" << _y0 << endl;
        cerr << "y0 ok" << endl;
        //_trajy.push_back(_y0);
    }
    else if (property=="Z0"){
        input >> _z0;
        coutf << "Z0" << _z0 << endl;
        cerr << "z0 ok" << endl;
        //_trajz.push_back(_z0);
    }
    else if (property=="METROPOLIS_STEP_0"){
        input >> _metropolis_step_0;
        coutf << "metropolis_step_0" << _metropolis_step_0 << endl;
        cerr << "step 0 ok" << endl;
    }
    else if (property=="X0_NEW"){
        input >> _x0_new;
        coutf << "X0_NEW" << _x0_new << endl;
        cerr << "x0 new ok" << endl;
        //_trajx.push_back(_x0);
    }
    else if (property=="Y0_NEW"){
        input >> _y0_new;
        coutf << "Y0_NEW" << _y0_new << endl;
        cerr << "y0 new ok" << endl;
        //_trajy.push_back(_y0);
    }
    else if (property=="Z0_NEW"){
        input >> _z0_new;
        coutf << "Z0_NEW" << _z0_new << endl;
        cerr << "z0 new ok" << endl;
        //_trajz.push_back(_z0);
    }
    else if (property=="METROPOLIS_STEP_1"){
        input >> _metropolis_step_1;
        coutf << "metropolis_step_1" << _metropolis_step_1 << endl;
        cerr << "step 1 ok" << endl;
    }
    else if (property=="SIGMA_0"){
        input >> _sigma_0;
        coutf << "sigma_0" << _sigma_0 << endl;
        cerr << "sigma 0 ok" << endl;
    }
    else if (property=="SIGMA_1"){
        input >> _sigma_1;
        coutf << "sigma_1" << _sigma_1 << endl;
        cerr << "sigma 1 ok" << endl;
    }
    else if( property == "NBLOCKS" ){
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << endl;
      cerr << "nblocks ok" << endl;
    } else if( property == "NSTEPS" ){
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << endl;
      cerr << "nsteps ok" << endl;
    } else if( property == "NITER" ){
      input >> _n_iter;
      coutf << "NITER= " << _n_iter << endl;
      cerr << "niter ok" << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}

void Metropolis :: initialize_properties(int state){
    string property;
    ifstream input("../INPUT/properties.dat");
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "r_medio" ){
                ofstream coutp("../OUTPUT/r.dat_" + to_string(state));
                coutp << "#     BLOCK:  ACTUAL_r:     r_AVE:      ERROR:" << endl;
                coutp.close();
                _measure_r = true ;
            }
            else if( property == "print_traj" ){
                ofstream coutj("../OUTPUT/traj.dat_" + to_string(state));
                coutj << "#         X:          Y:          Z:" << endl;
                coutj.close();
                _print_r = true;
                break;
            }
            else if( property == "ENDPROPERTIES" ){
                ofstream coutf;
                coutf.open("../OUTPUT/output.dat",ios::app);
                coutf << "Reading properties completed!" << endl;
                coutf.close();
                break;
            } else cerr << "PROBLEM: unknown property" << endl;
        }
    input.close();
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl;
}

int Metropolis::get_nbl() {
    return _nblocks;
}

int Metropolis::get_nsteps() {
    return _nsteps;
}

int Metropolis::get_nbl_eq() {
    return _n_eq_blocks;
}

int Metropolis::get_nsteps_eq() {
    return _n_eq_steps;
}

int Metropolis::get_nattempts() {
    return _nattempts;
}

int Metropolis::get_naccepted() {
    return _naccepted;
}

double Metropolis::mod_psi_0(double x, double y, double z) {
    double r_squared = x*x + y*y + z*z;
    double exponent = -2.0 * sqrt(r_squared) ;
    return (1.0 / (M_PI)) * exp(exponent);
}



//double Metropolis::mod_psi_1(double x, double y, double z, double a0) {
  //  double r_squared = x*x + y*y + z*z;
   // double r = sqrt(r_squared);
   // double exponent = -r / (2.0 * a0);
    //return (1.0 / (32.0 * M_PI * a0 * a0 * a0 * a0 * a0)) * r_squared * exp(exponent);
//}

//double Metropolis :: mod_psi_1(double x, double y, double z) {
    // Costanti
   //double coeff = 1.0 / 64.0;
    //double sqrt2_over_pi = sqrt(2.0 / M_PI);
    //double norm_squared = x * x + y * y + z * z;

    // Calcolo del modulo quadro
   // double result = coeff * sqrt2_over_pi * norm_squared * exp(-sqrt(norm_squared));
   // return result;
//}

double Metropolis :: mod_psi_1(double x, double y, double z) {
    // Costanti
    double coeff = 1.0 / (32.0* M_PI);
    double r_squared = x * x + y * y + z * z;

    // Calcolo del modulo quadro
    double result = coeff  * exp(-sqrt(r_squared)) * z * z ;
    return result;
}

void Metropolis::newMetropolis_algorithm(int state) {
    vec position;
    position.resize(_ndim);
    position(0) = _x;
    position(1) = _y;
    position(2) = _z;
    double x,y,z,xnew,ynew,znew,alfa,extraction;
    //for(int i=0; i< _n_iter; i++){
        // x=x_n
        x = position(0);
        y = position(1);
        z = position(2);
        _trajx = x;
        _trajy = y;
        _trajz = z ;
        //_particle.setposition(0, x);
        //_particle.setposition(1, y);
        //_particle.setposition(2, z);
        _nattempts++;

        // generate x' from T(x'|x_n) and evaluate A(x'|x_n) = alfa
        if (state == 0) {
            xnew = _rnd.Rannyu(x-1.0*_metropolis_step_0,x+1.0*_metropolis_step_0);
            ynew = _rnd.Rannyu(y-1.0*_metropolis_step_0,y+1.0*_metropolis_step_0);
            znew = _rnd.Rannyu(z-1.0*_metropolis_step_0,z+1.0*_metropolis_step_0);
            alfa = std::min(1.,mod_psi_0(xnew,ynew,znew)/mod_psi_0(x,y,z));
        }
        
        else if (state == 1){
            xnew = _rnd.Rannyu(x-1.0*_metropolis_step_1,x+1.0*_metropolis_step_1);
            ynew = _rnd.Rannyu(y-1.0*_metropolis_step_1,y+1.0*_metropolis_step_1);
            znew = _rnd.Rannyu(z-1.0*_metropolis_step_1,z+1.0*_metropolis_step_1);
            alfa = std::min(1.,mod_psi_1(xnew,ynew,znew)/mod_psi_1(x,y,z));
        }

        // accept the move with prob alfa
        extraction = _rnd.Rannyu(0,1);
        if(extraction <= alfa){    //accept
            position(0) = xnew;
            position(1) = ynew;
            position(2) = znew;
            _naccepted++; 
            this->set_xyz(xnew,ynew,znew);
        }
        _block_av += sqrt(position(0)*position(0)+position(1)*position(1)+position(2)*position(2));  //aggiungo r al blocco
        // repeat
    //}
    //_block_av = _block_av/_n_iter;
    //return;   
}

void Metropolis::Metropolis_algorithm(double x0, double y0, double z0, int state) {
    vec position;
    position.resize(_ndim);
    position(0) = x0;
    position(1) = y0;
    position(2) = z0;
    double x,y,z,xnew,ynew,znew,alfa,extraction;
    //for(int i=0; i< _n_iter; i++){
        // x=x_n
        x = position(0);
        y = position(1);
        z = position(2);
        _trajx = x;
        _trajy = y;
        _trajz = z ;
        //_particle.setposition(0, x);
        //_particle.setposition(1, y);
        //_particle.setposition(2, z);
        _nattempts++;

        // generate x' from T(x'|x_n) and evaluate A(x'|x_n) = alfa
        if (state == 0) {
            xnew = _rnd.Rannyu(x-1.0*_metropolis_step_0,x+1.0*_metropolis_step_0);
            ynew = _rnd.Rannyu(y-1.0*_metropolis_step_0,y+1.0*_metropolis_step_0);
            znew = _rnd.Rannyu(z-1.0*_metropolis_step_0,z+1.0*_metropolis_step_0);
            alfa = std::min(1.,mod_psi_0(xnew,ynew,znew)/mod_psi_0(x,y,z));
        }
        
        else if (state == 1){
            xnew = _rnd.Rannyu(x-1.0*_metropolis_step_1,x+1.0*_metropolis_step_1);
            ynew = _rnd.Rannyu(y-1.0*_metropolis_step_1,y+1.0*_metropolis_step_1);
            znew = _rnd.Rannyu(z-1.0*_metropolis_step_1,z+1.0*_metropolis_step_1);
            alfa = std::min(1.,mod_psi_1(xnew,ynew,znew)/mod_psi_1(x,y,z));
        }

        // accept the move with prob alfa
        extraction = _rnd.Rannyu(0,1);
        if(extraction <= alfa){    //accept
            position(0) = xnew;
            position(1) = ynew;
            position(2) = znew;
            _naccepted++; 
            
        }
        _block_av += sqrt(position(0)*position(0)+position(1)*position(1)+position(2)*position(2));  //aggiungo r al blocco
        // repeat
    //}
    //_block_av = _block_av/_n_iter;
    //return;   
}

void Metropolis::Metropolis_algorithm_Gauss(int state) {
    vec position;
    position.resize(_ndim);
    position(0) = _x;
    position(1) = _y;
    position(2) = _z;
    double x,y,z,xnew,ynew,znew,alfa,extraction;
    //for(int i=0; i< _n_iter; i++){
        // x=x_n
        x = position(0);
        y = position(1);
        z = position(2);
        _trajx = x;
        _trajy = y;
        _trajz = z ;
        //_particle.setposition(0, x);
        //_particle.setposition(1, y);
        //_particle.setposition(2, z);
        _nattempts++;

        // generate x' from T(x'|x_n) and evaluate A(x'|x_n) = alfa
        if (state == 0) {
            xnew = _rnd.Gauss(x,_sigma_0);
            ynew = _rnd.Gauss(y,_sigma_0);
            znew = _rnd.Gauss(z,_sigma_0);
            alfa = std::min(1.,mod_psi_0(xnew,ynew,znew)/mod_psi_0(x,y,z));
        }
        
        else if (state == 1){
            xnew = _rnd.Gauss(x,_sigma_1);
            ynew = _rnd.Gauss(y,_sigma_1);
            znew = _rnd.Gauss(z,_sigma_1);
            alfa = std::min(1.,mod_psi_1(xnew,ynew,znew)/mod_psi_1(x,y,z));
        }

        // accept the move with prob alfa
        extraction = _rnd.Rannyu(0,1);
        if(extraction <= alfa){    //accept
            position(0) = xnew;
            position(1) = ynew;
            position(2) = znew;
            _naccepted++; 
            this->set_xyz(xnew,ynew,znew);
        }
        _block_av += sqrt(position(0)*position(0)+position(1)*position(1)+position(2)*position(2));  //aggiungo r al blocco
        // repeat
    //}
    //_block_av = _block_av/_n_iter;
    //return;   
}


void Metropolis::r_averages(int blk, int state){
    
    ofstream coutf;
    //_average     = _block_av /double(_nsteps) /double(_n_iter);
    _average     = _block_av /double(_nsteps);
    _global_av  += _average;
    _global_av2 += _average * _average; // % -> element-wise multiplication
    
    if (_measure_r) {
        coutf.open("../OUTPUT/r.dat_" + to_string(state),ios::app);
        coutf << setw(15) << blk 
            << setw(15) << _average
            << setw(15) << _global_av/double(blk)
            << setw(15) << this->error(_global_av, _global_av2, blk) << endl;
        coutf.close();
    }

    return;
}

void Metropolis :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av = 0 ;
  return;
}

double Metropolis :: get_initial_x(){
    return _x0;
}

double Metropolis :: get_initial_y(){
    return _y0;
}

double Metropolis :: get_initial_z(){
    return _z0;
}

double Metropolis :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

void Metropolis :: set_nattempts(){
    _nattempts = 0;
}

void Metropolis :: set_naccepted(){
    _naccepted = 0;
}

void Metropolis :: set_r0(){
    _x0 = _x0_new;
    _y0 = _y0_new;
    _z0 = _z0_new;
}

void Metropolis :: set_xyz(double x,double y, double z){
    _x = x;
    _y = y;
    _z = z;
}



void Metropolis :: print_traj(int state){
    ofstream coutj;
    if (_print_r) {
        coutj.open("../OUTPUT/traj.dat_" + to_string(state),ios::app);
        coutj << setw(15) << _trajx 
            << setw(15) << _trajy
            << setw(15) << _trajz << endl;
        coutj.close();
    }

    return;
}

