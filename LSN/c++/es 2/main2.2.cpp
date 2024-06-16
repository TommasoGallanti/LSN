#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "random.h"
#include "blockfunctions.h"
#include "RWalk.h"

using namespace std;

int main(int argc, char *argv[]) {


    Random rnd;
    rnd.StartGen();


    vector<double> rndm(300);
    vector<double> rndm_c(300);
    int N{100}; //  blocks
    int M{100}; //  walks for block
    int steps{100};

    Walk value(steps);
    Walk value_c{steps};

    vector<double> err_prog(N), r2(N), r_blk(N), r_blk2(N);
    vector<double> err_prog_c(N), r2_c(N), r_blk_c(N), r_blk2_c(N);


        //discreto//

    for(int l{}; l < N ; ++l) {

        r2.assign(M,0.);    //start each walk per block from zero

        for(int k{}; k < M; ++k){ //create a block

            for(size_t i{}; i < rndm.size(); ++i)
                rndm[i] = rnd.Rannyu();
        
            value.Random_Walk(rndm);

            for(size_t j{}; j < r2.size(); ++j){
                    r2[j] += value.GetR2(j)/double(M);

        
            }
        }

        for(int i{}; i < N; ++i) {
            r_blk[i] += sqrt(r2[i])/double(N);
            r_blk2[i] += r2[i]/double(N);
        }

        for(int i{}; i < N; ++i) 
            err_prog[i] = error(r_blk,r_blk2,i);

    }

    ofstream out_discrete("Discrete_RW.dat");

    for(int i{}; i < N; ++i) 
        out_discrete<<r_blk[i] <<" "<<err_prog[i]<<endl;


    
        //continuo//
    
    for(int l{}; l < N ; ++l) {

        r2_c.assign(M,0.);    //start each walk per block from zero

        for(int k{}; k < M; ++k){ //create a block

            for(size_t i{}; i < rndm_c.size(); ++i)
                rndm_c[i] = rnd.Rannyu();
        
            value_c.Random_Cont(rndm_c);

            for(size_t j{}; j < r2_c.size(); ++j){
                    r2_c[j] += value_c.GetR2(j)/double(M);

        
            }
        }

        for(int i{}; i < N; ++i) {
            r_blk_c[i] += sqrt(r2_c[i])/double(N);
            r_blk2_c[i] += r2_c[i]/double(N);
        }

        for(int i{}; i < N; ++i) 
            err_prog_c[i] = error(r_blk_c,r_blk2_c,i);

    }

    ofstream out_continue("Continue_RW.dat");

    for(int i{}; i < N; ++i) {
        out_continue<<r_blk_c[i] <<" "<<err_prog_c[i]<<endl;

    }

    return 0;
}