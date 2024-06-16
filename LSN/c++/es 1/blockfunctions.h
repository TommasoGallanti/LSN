#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

double error(vector<double> mean, vector<double> mean2, int n) { 
    return (n == 0 ? 0 : sqrt((mean2[n] - (mean[n]*mean[n]))/double(n)));
}

void Block_statistic( vector<double> r , int n_blocks, int n_elements, string outFile) {

    int N = n_blocks;
    int M = n_elements;
    int L = M/N;

    vector<double> mean(N), mean2(N), sum_prog(N), sum2_prog(N),err_prog(N);

    for(auto i=0; i < N ; ++i ) {
      
        double sum{};

        for(auto j{0}; j < L ; ++j ) {
            int k = j + i*L;
            sum += r[k]; 
        }
        mean[i] = sum/double(L);
        mean2[i] = mean[i]*mean[i];
    }

   ofstream out(outFile);
   
    for(auto i=0; i < N ; ++i ) {
        for(auto j=0; j < i+1 ; ++j ) {
            sum_prog[i] += mean[j];
            sum2_prog[i] += mean2[j];
        }
        sum_prog[i]/=(i+1);
        sum2_prog[i]/=(i+1);
        err_prog[i] = error(sum_prog, sum2_prog,i);

        out << sum_prog[i] <<" "<< err_prog[i]<<endl;
   }

   


}