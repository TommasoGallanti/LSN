
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "random.h"
#include "blockfunctions.h"

using namespace std;

double chi(double ni, int n, int M) {
   return (n == 0 ? 0 : (pow(ni - n / M, 2) / (n / M)));
}

int main(int argc, char *argv[])
{

   Random rnd;
   rnd.StartGen();

   // es. 0.1.1.1 mean r

   int M = 100000, N = 100;

   vector<double> r; 
   for (auto i = 0; i < M; ++i)   
      r.push_back(rnd.Rannyu(0., 1.));

   Block_statistic(r,N,M,"ex_1.1.1.dat");

   // es. 0.1.1.2 

   r.resize(0);
   for (auto i = 0; i < M; ++i)
      r.push_back( pow(rnd.Rannyu() -0.5 , 2) ); 

   Block_statistic(r,N,M,"ex_1.1.2.dat");


   // es.0.1.1.3 chi squared

   M=100;
   int n= 10000;
  
   ofstream outChi2("ex_1.1.3.dat");

   double chi2_value{};

   vector<double> chi_progres(M);

   for(int k=0; k<100; ++k){
      chi2_value = 0;
      for(int i=0; i < M; ++i) {
         int count=0;
         for(int j=0; j < n; ++j) {
            double rndmNum=rnd.Rannyu();
            if( (i)/double(M)<= rndmNum && rndmNum<=(i+1)/double(M) ){ //verifico se rndmNum appartiene all'intervallo i-esimo
               count++;
            }
         } 
         chi_progres[i]=chi(count,n,M);
         chi2_value+= chi_progres[i];
      }
      outChi2<<chi2_value<<endl;
   }
   return 0;
}