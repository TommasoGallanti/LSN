#include <iostream>
#include "varMC.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main (int argc, char *argv[]){

    Var_MC Var;
    Var.Initialize();
    Var.initialize_properties();
    
    //create temperatures vector
    const int size = Var.get_n_temps();;
    double temperatures[size];

    double up = Var.get_temp_up();
    double low = Var.get_temp_low();

    double step = (up - low) / (size);

    for (int i = 0; i < size+1; ++i) {
        temperatures[i] = up - i * step;
    }
    temperatures[size] = low;


    cout << "t final = " << temperatures[size] << endl;
    
    cout << "le temperature per i SA sono: " << endl;
    
    for (int i = 0; i < size+1; ++i) {
        cout << i << "  " << temperatures[i] << endl;
    }

    //START SA ALGORITHM 
    for(int l=0; l < size ; l++){
      Var.set_temp(temperatures[l]);
      Var.SimulatedAnnealing();
    }

    cout << "accepted = " << Var.get_naccepted() << endl;
    cout << "attemps  = " << Var.get_nattempts() << endl;
    cout << "ratio    = " << double(Var.get_naccepted())/double(Var.get_nattempts()) << endl;

    
    Var.set_best_parameters("../OUTPUT/Results.dat");      //imposta mu e sigma come le migliori trovate durante l'ottimizzazione
    cout << "best mu    = " << Var.get_mu() << endl;
    cout << "best sigma = " << Var.get_sigma() << endl;
   
   
    //DATA BLOCKING
    ofstream outx;
    outx.open("../OUTPUT/traj.dat");
    for (int i=0; i<Var.get_nbl(); i++){
        for (int j=0; j<Var.get_nsteps(); j++){
            Var.Metropolis_algorithm(Var.get_mu(),Var.get_sigma());
            outx << Var.get_x() << endl;
            Var.measure();
        }
        Var.averages(i+1);
        cout << "blocco" << i+1 << endl;
        Var.block_reset(i+1);
        
    }
    outx.close();



  cout<<"programma eseguito"<<endl;
  return 0;
}



