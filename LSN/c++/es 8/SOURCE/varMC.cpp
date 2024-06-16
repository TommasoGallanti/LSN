#include <cmath>
#include <cstdlib>
#include <string>
#include "varMC.h"
#include <armadillo>

using namespace std;
using namespace arma;


double Var_MC:: Gaussiana(double x, double mu, double sigma) {
    return exp(-pow( (x-mu)/sigma, 2.) /2.);
}

double Var_MC::Kinetic(double x, double sigma, double mu) {
    
    double a = (x-mu)*(x-mu)/(sigma*sigma);
    double b = (x+mu)*(x+mu)/(sigma*sigma);
    return 0.5/(sigma*sigma)*( 1. -( a*Gaussiana(x,mu,sigma) + b*Gaussiana(x,-mu,sigma) )/Psi(x,sigma,mu) );
}

double Var_MC::Potenzial(double x) {
    return pow(x,4) -2.5*pow(x,2);
}
double Var_MC::Hamiltonian(double x, double sigma,double mu) {
    return Potenzial(x) + Kinetic(x,sigma, mu);
}

double Var_MC:: Psi(double x,double sigma, double mu) {
 return Gaussiana(x,mu,sigma) + Gaussiana(x,-mu,sigma);
}
double Var_MC::mod_psi(double x,double sigma,double mu){
    return Psi(x,sigma,mu)*Psi(x,sigma,mu) ;
}

void Var_MC::Metropolis_algorithm(double mu, double sigma) {
    
    double x = _x;
    double xnew,alfa,extraction;
    
    // x=x_n
    //_trajx = x;
    _nattempts++;

        // generate x' from T(x'|x_n) and evaluate A(x'|x_n) = alfa
        
    //xnew = _rnd.Rannyu(x-1.0*_metropolis_step,x+1.0*_metropolis_step);
    xnew =(x + _metropolis_step*(_rnd.Rannyu() -0.5));
    alfa = std::min(1.,mod_psi(xnew,sigma,mu)/mod_psi(x,sigma,mu));
        
        

        // accept the move with prob alfa
    extraction = _rnd.Rannyu(0,1);
    if(extraction <= alfa){    //accept
        _x = xnew;
        //_mu = munew;
        //_sigma = sigmanew;

        _naccepted++; 
        }


    //_block_av += this->Hamiltonian();  //aggiungo x al blocco
    // repeat
    
    //_block_av = _block_av/_n_iter;
    return;   
}

void Var_MC :: SimulatedAnnealing(){
    double H_old = 0;
    double H_new = 0;
    double mu_new, sigma_new;
    double x = _x;
    int wd = 20;


    ofstream coutp;
    coutp.open("../OUTPUT/Results.dat",ios::app);
    
    for (int k{}; k < _temp_rep; k++){
        H_old = 0;
        H_new = 0;
        _x = x;
        
        for (int i=0; i<_n_iter; i++){
            this->Metropolis_algorithm(_mu,_sigma); 
            H_old += Hamiltonian(_x,_sigma,_mu);
        }
        H_old /= _n_iter;
        mu_new =fabs( _mu + _delta_opt*(_rnd.Rannyu() -0.5) );     //estrazione di mu e sigma
        sigma_new =fabs( _sigma + _delta_opt*(_rnd.Rannyu() - 0.5) ); 
        _x = x; //Reset of the starting position

        for(int i=0 ; i< _n_iter; ++i){
            this->Metropolis_algorithm(mu_new,sigma_new);
            H_new += Hamiltonian(_x,sigma_new,mu_new);
        }
        H_new /= _n_iter;
                
                

        double p = exp( 1./_temp*(H_old-H_new)) ;
        double a = _rnd.Rannyu();
                
        if(p >= a ){
            _mu = mu_new;
            _sigma = sigma_new;
            _SA_accepted++;
            
        }
        _SA_attempted++;
    }
    coutp << _temp<< setw(wd) << _mu << setw(wd)<< _sigma << setw(wd) << H_new <<endl;
    coutp.close();
}


void Var_MC :: initialize_properties(){
    _nprop=0;
    int index = 0;
    _measure_penergy  = false; //Defining which properties will be measured
    _measure_kenergy  = false;
    _measure_tenergy  = false;
    string property;
    ifstream input("../INPUT/properties.dat");
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "mu_sigma_H" ){
                ofstream coutp("../OUTPUT/Results.dat");
                coutp << "#     Temp:      mu:     sigma:      H:" << endl;
                coutp.close();
                cerr << "prop mu sigma H" << endl;
            }
            else if( property == "Potential" ){
                ofstream coutj("../OUTPUT/Potential.dat");
                coutj << "#     BLOCK:    ACTUAL_PE:     PE_AVE:      ERROR: " << endl;
                coutj.close();
                _nprop++;
                _measure_penergy=true;
                _index_penergy = index;
                index++;
                
            }
            else if( property == "Kinetic" ){
                ofstream coutj("../OUTPUT/Kinetic.dat");
                coutj << "#     BLOCK:    ACTUAL_PE:     PE_AVE:      ERROR: " << endl;
                coutj.close();
                _nprop++;
                _measure_kenergy=true;
                _index_kenergy=index;
                index++;
            }
            else if( property == "Total" ){
                ofstream coutj("../OUTPUT/H.dat");
                coutj << "#     BLOCK:    ACTUAL_PE:     PE_AVE:      ERROR: " << endl;
                coutj.close();
                _nprop++;
                _measure_tenergy=true;
                _index_tenergy=index;
                index++;
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

  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
}


void Var_MC :: Initialize(){ // Initialize the Var_MC object according to the content of the input files in the ../INPUT/ directory

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
        input >> _x;
        coutf << "X0 = " << _x << endl;
        cerr << "x0 ok" << endl;
        //_trajx.push_back(_x0);
    }
    else if (property=="METROPOLIS_STEP"){
        input >> _metropolis_step;
        coutf << "metropolis_step = " << _metropolis_step << endl;
        cerr << "metro step ok" << endl;
    }

    else if (property=="DELTA_OPT"){
        input >> _delta_opt;
        coutf << "delta_opt = " << _delta_opt << endl;
        cerr << "delta opt ok" << endl;
    }
    else if (property=="SIGMA_0"){
        input >> _sigma;
        coutf << "sigma_0 = " << _sigma << endl;
        cerr << "sigma 0 ok" << endl;
    }

    else if (property=="MU_0"){
        input >> _mu;
        coutf << "mu_0 = " << _mu << endl;
        cerr << "mu 0 ok" << endl;
    }
    
    else if (property=="TEMP"){
        //input >> _temp;
        input >> _t_low;
        input >> _t_up;
        input >> _n_temps;
        cerr << "temp ok" << endl;
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
    } else if( property == "TEMP_REP" ){
      input >> _temp_rep;
      coutf << "temp rep= " << _temp_rep << endl;
      cerr << "temp rep ok" << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  coutf << "Var_MC initialized!" << endl;
  coutf.close();
  return;
    
}



void Var_MC :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av.zeros();
  
  return;
}

int Var_MC :: get_nbl(){
  return _nblocks;
}

int Var_MC :: get_nsteps(){
  return _nsteps;
}


void Var_MC :: set_temp(double temp){
  _temp = temp;
}
double Var_MC :: get_temp(){
  return _temp;
}



int Var_MC::get_nattempts() {
    return _nattempts;
}

int Var_MC::get_naccepted() {
    return _naccepted;
}

int Var_MC::get_n_temps() {
    return _n_temps;
}

double Var_MC::get_temp_up(){
    return _t_up;
}

double Var_MC::get_temp_low(){
    return _t_low;
}

double Var_MC::get_mu(){
    return _mu;
}

double Var_MC::get_sigma(){
    return _sigma;
}

void Var_MC::set_mu(double mu){
    _mu = mu;
}

void Var_MC::set_sigma(double sigma){
    _sigma = sigma;
}


void Var_MC::set_best_parameters(string filename){
    int k = _n_temps;
    double temp,mu,sigma,H,H_min;
    ifstream in;
    in.open(filename);
    
    // Skip the header line
    std::string header;
    std::getline(in, header);

    in>>temp>>mu>>sigma>>H_min;
    
    for(int i{1}; i < k+1; ++i){
        in>>temp>>mu>>sigma>>H;
        if(H < H_min) {
            _mu = mu;
            _sigma = sigma;
            H_min = H;
        }
    } 
    in.close();  
}

void Var_MC::averages(int blk){
    
    ofstream coutf;
    double average, sum_average, sum_ave2;
    _average     = _block_av /double(_nsteps);
    _global_av  += _average;
    _global_av2 += _average % _average; // % -> element-wise multiplication
    
    
    
    if (_measure_tenergy) {
        coutf.open("../OUTPUT/H.dat" ,ios::app);
        average  = _average(_index_tenergy);
        sum_average = _global_av(_index_tenergy);
        sum_ave2 = _global_av2(_index_tenergy);
        coutf << setw(15) << blk 
            << setw(15) << average
            << setw(15) << sum_average/double(blk)
            << setw(15) << this->error(sum_average, sum_ave2, blk) << endl;
        coutf.close();
    }

    if (_measure_penergy) {
        coutf.open("../OUTPUT/Potential.dat" ,ios::app);
        average  = _average(_index_penergy);
        sum_average = _global_av(_index_penergy);
        sum_ave2 = _global_av2(_index_penergy);
        coutf << setw(15) << blk 
            << setw(15) << average
            << setw(15) << sum_average/double(blk)
            << setw(15) << this->error(sum_average, sum_ave2, blk) << endl;
        coutf.close();
    }

    if (_measure_kenergy) {
        coutf.open("../OUTPUT/Kinetic.dat" ,ios::app);
        average  = _average(_index_kenergy);
        sum_average = _global_av(_index_kenergy);
        sum_ave2 = _global_av2(_index_kenergy);
        coutf << setw(15) << blk 
            << setw(15) << average
            << setw(15) << sum_average/double(blk)
            << setw(15) << this->error(sum_average, sum_ave2, blk) << endl;
        coutf.close();
    }

    return;
}



void Var_MC :: measure(){ // Measure properties
  
  // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
  double penergy_temp=0.0; // temporary accumulator for potential energy
  double kenergy_temp=0.0; // temporary accumulator for kinetic energy
  double tenergy_temp=0.0;

  
  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    //penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
    //double x = sqrt(pow(_particle.getposition(0,true),2)+pow(_particle.getposition(1,true),2)+pow(_particle.getposition(2,true),2));
    penergy_temp = this->Potenzial(_x);
    _measurement(_index_penergy) = penergy_temp;
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    //for (int i=0; i<_npart; i++) kenergy_temp += 0.5 * dot( _particle(i).getvelocity() , _particle(i).getvelocity() ); 
    kenergy_temp = Kinetic(_x,_sigma,_mu);
    _measurement(_index_kenergy) = kenergy_temp;
  }
  // TOTAL ENERGY (kinetic+potential) //////////////////////////////////////////
  if (_measure_tenergy){
     _measurement(_index_tenergy) = kenergy_temp + penergy_temp;
  }
  

  _block_av += _measurement; //Update block accumulators
  

  return;
}

double Var_MC :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

double Var_MC::get_x(){
    return _x;
}