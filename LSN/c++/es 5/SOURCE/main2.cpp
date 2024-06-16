#include <iostream>
#include "metropolis.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main (int argc, char *argv[]){

    //double a0 = 1.;    // setto raggio di Bohr come unità
    double x0,y0,z0;
    Metropolis MTR;
    int state=1;
    MTR.Initialize(state);
    cout << MTR.get_nattempts() << endl;
    x0 = MTR.get_initial_x();
    y0 = MTR.get_initial_y();
    z0 = MTR.get_initial_z();

    cout << "calcolo valor medio su stato 100" << endl;
    for (int i=0; i<MTR.get_nbl(); i++){
        for (int j=0; j<MTR.get_nsteps(); j++){
            MTR.Metropolis_algorithm(x0,y0,z0,state);
        }
        MTR.r_averages(i+1,state);
        cout << "blocco" << i+1 << endl;
        MTR.block_reset(i+1);
    }
    double n_attemps = MTR.get_nattempts();
    double n_accepted = MTR.get_naccepted();
    cout << "n_attemps = "<< MTR.get_nattempts() << endl;
    cout << "n_accepted = " << MTR.get_naccepted() << endl;
    cout << "acceptance ratio = " << n_accepted/n_attemps << endl;


   


return 0;

}