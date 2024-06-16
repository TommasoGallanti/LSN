#include <iostream>
#include "metropolis.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main (int argc, char *argv[]){

    //double a0 = 1.;    // setto raggio di Bohr come unità
    double x0,y0,z0;
    Metropolis MTR;
    int state=0;
    double n_attemps,n_accepted;
    MTR.Initialize(state);
    MTR.initialize_properties(state);
    MTR.set_nattempts();
    cout << MTR.get_nattempts() << endl;
    x0 = MTR.get_initial_x();
    y0 = MTR.get_initial_y();
    z0 = MTR.get_initial_z();
    MTR.set_xyz(x0,y0,z0);

    cout << "calcolo valor medio su stato 100" << endl;
    for (int i=0; i<MTR.get_nbl(); i++){
        for (int j=0; j<MTR.get_nsteps(); j++){
            
            //MTR.newMetropolis_algorithm(state);
            //MTR.Metropolis_algorithm(x0,y0,z0,state);
            MTR.Metropolis_algorithm_Gauss(state);
            if(j%100 == 0){
            MTR.print_traj(state);
        }
        }
        MTR.r_averages(i+1,state);
        cout << "blocco" << i+1 << endl;
        MTR.block_reset(i+1);
    }
    n_attemps = MTR.get_nattempts();
    n_accepted = MTR.get_naccepted();
    cout << "n_attemps = "<< MTR.get_nattempts() << endl;
    cout << "n_accepted = " << MTR.get_naccepted() << endl;
    cout << "acceptance ratio = " << n_accepted/n_attemps << endl;

    
    Metropolis MTR2;
    state=1;
    MTR2.Initialize(state);
    MTR2.initialize_properties(state);
    MTR2.set_nattempts();
    cout << MTR2.get_nattempts() << endl;
    x0 = MTR2.get_initial_x();
    y0 = MTR2.get_initial_y();
    z0 = MTR2.get_initial_z();
    MTR2.set_xyz(x0,y0,z0);

    cout << "calcolo valor medio su stato 210" << endl;
    for (int i=0; i<MTR2.get_nbl(); i++){
        for (int j=0; j<MTR2.get_nsteps(); j++){
            //MTR2.newMetropolis_algorithm(state);
            //MTR2.Metropolis_algorithm(x0,y0,z0,state);
            MTR2.Metropolis_algorithm_Gauss(state);
            if(j%100 == 0){
            MTR2.print_traj(state);
        }
        }
        MTR2.r_averages(i+1,state);
        cout << "blocco" << i+1 << endl;
        MTR2.block_reset(i+1);
        
    }
    n_attemps = MTR2.get_nattempts();
    n_accepted = MTR2.get_naccepted();
    cout << "n_attemps = "<< MTR2.get_nattempts() << endl;
    cout << "n_accepted = " << MTR2.get_naccepted() << endl;
    cout << "acceptance ratio = " << n_accepted/n_attemps << endl;


    


return 0;

}