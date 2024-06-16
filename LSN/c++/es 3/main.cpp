#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "random.h"
#include "blockfunctions.h"
#include "funzioni_statistiche.h"

using namespace std;

double S(double t, double S0, double mu, double z, double sigma) {
    return S0 * exp((mu - sigma * sigma / 2.0) * t + sigma * z * sqrt(t));
}


int main(int argc, char *argv[]) {


    Random rnd;
    rnd.StartGen();

    
    double S0=100;
    double T=1;
    double K=100;
    double mu=.1;
    double sigma=.25;

    int N{100}; //  blocks
    int M{10000}; // iteration 

    double rndn;

    


    //direct calculation
    double S_dir{};
    vector<double> C_dir;
    vector<double> P_dir;
    
    for(int i{}; i < M; ++i) {     //only one step 
        rndn = rnd.Gauss(0.,1.);
        S_dir = S(T,S0,mu,rndn,sigma);
        C_dir.push_back(exp(-mu*T)*max(0.,S_dir -K));
        P_dir.push_back(exp(-mu*T)*max(0.,K - S_dir));
    }


    Block_statistic(C_dir,N,M,"Call_direct.dat");
    Block_statistic(P_dir,N,M,"Put_direct.dat");

    double Call_mean = CalcolaMedia(C_dir);
    double Put_mean = CalcolaMedia(P_dir);


    cout<< "Direct method:    "<<"Call_mean "<<Call_mean<<" "<<"Put_mean "<<Put_mean<<endl;


    //discrete calculation
    int steps{100};
    double deltat= T/double(steps);
    vector<double> C_discr;
    vector<double> P_discr;

    for(int i{}; i < M; ++i) {     
        
        double S_iter{S0};
        double S_discr;
        
        
        for(int j{}; j<steps; j++) {   //100 steps 
            rndn = rnd.Gauss(0.,1.);
            S_discr = S(deltat,S_iter,mu,rndn,sigma);
            S_iter=S_discr;
        } 

        C_discr.push_back(exp(-mu*T)*max(0.,S_discr -K));
        P_discr.push_back(exp(-mu*T)*max(0.,K - S_discr));
    }

    Block_statistic(C_discr,N,M,"Call_discr.dat");
    Block_statistic(P_discr,N,M,"Put_discr.dat");

    cout<< "Discrete method:  "<<"Call_mean "<<Call_mean<<" "<<"Put_mean "<<Put_mean<<endl;

    return 0;


}