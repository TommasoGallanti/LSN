#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "random.h"
#include "blockfunctions.h"

using namespace std;

#include <cmath>


double func(double x) {
    return M_PI / 2.0 * cos(M_PI * x / 2.0);
}


double p(double x) {
    return 2 * (1-x);
}


double g(double x) {
    return func(x) / p(x);
}


double inv(double y) {
    return 1 - sqrt(1 - y);
}


int main(int argc, char *argv[])
{

    Random rnd;
    rnd.StartGen();
     
    int M = 1000000;
    int N = 100;

    vector<double> media(N), media2(N), sum_prog(N), sum2_prog(N), err_prog(N);
    vector<double> mediaT(N), media2T(N), sumT_prog(N), sum2T_prog(N), errT_prog(N);

    for(auto i{0}; i<N; i++) {
        
        double sum{};
        double sumT{};

        for(auto j{0}; j<M/N; j++) {
            sum+= func(rnd.Rannyu());
            sumT+=g(inv(rnd.Rannyu()));
        }

        media[i] = sum/double(M/N);
        mediaT[i] = sumT/double(M/N);
        media2[i] = media[i]*media[i];
        media2T[i] = mediaT[i]*mediaT[i];
    }
        
    ofstream out1("Unif.dat");
    ofstream out2("IS.dat");

    for(auto i{0}; i < N; ++i) {
        for(auto j{0}; j < i+1 ; ++j ) {
            sum_prog[i] += media[j];
            sum2_prog[i] += media2[j];

            sumT_prog[i] += mediaT[j];
            sum2T_prog[i] += media2T[j];
        }
        sum_prog[i]/=(i+1);
        out1 << sum_prog[i] << " ";

        sum2_prog[i]/=(i+1);

        err_prog[i] = error(sum_prog,sum2_prog,i);
        out1 << err_prog[i]<<endl;

        sumT_prog[i]/=(i+1);
        out2 << sumT_prog[i] << " ";

        sum2T_prog[i]/=(i+1);

        errT_prog[i] = error(sumT_prog,sum2T_prog,i);
        out2 << errT_prog[i]<<endl;
    }

  
return 0;





}