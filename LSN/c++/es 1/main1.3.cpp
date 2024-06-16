#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "random.h"
#include "blockfunctions.h"

using namespace std;



int main(int argc, char *argv[])
{
   Random rnd;

   rnd.StartGen();

     

   vector<double> v(2);
   double L = 7; 
   double d = 15; 

   int M = 1000000; //number of throws
   int N = 100;  //number of blocks
   int K = M/N ; 


   vector<double> pi(N), pi2(N), error_prog(N), pi_prog(N), pi2_prog(N); 

   ofstream piOut("pi.dat");

   double theta{};

   //cout<<"inizia primo for"<<endl;

   for(int j=0; j < N; ++j){
      double hit{}, tot{};
      //cout<<j<<endl;
      //cout<<"primo for concat"<<endl;
      for(int i=0; i < K; ++i) {
         v[0] = rnd.Rannyu(0.,double(1000*d));   //initial point 
         //cout<<"v[0]="<<v[0]<<endl;
         theta = rnd.RanAngle();
         //cout<<theta<<endl;
         v[1] = v[0] + cos(theta)*L;    //projection whit costheta
         //cout<<i<<endl;
         
         if(int(v[0]/d) != int(v[1]/d))
            hit++;
         tot++;
      }

      double P= hit/double(tot);
      pi[j] = 2*L / (P*d);
      pi2[j] = pi[j]* pi[j];
   }

    //cout<<"inizia secondo for"<<endl;

   for(int i=0; i < N ; ++i ) {
      for(int j=0; j < i+1 ; ++j ) {

         pi_prog[i] += pi[j];
         pi2_prog[i] += pi2[j];
      }

      pi_prog[i]/=(i+1);
      pi2_prog[i]/=(i+1);
      error_prog[i] = error(pi_prog,pi2_prog,i);

      //cout<<i<<endl;
      piOut<<pi[i]<<" ";
      piOut<<pi_prog[i]<<" ";
      piOut<<error_prog[i]<<endl;
   }

   return 0;
}


