

#include <iostream>
#include "system.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main (int argc, char *argv[]){

    const int size = 15;
    double temperatures[size];

    double start = 0.5;
    double end = 2.0;

    double step = (end - start) / (size - 1);

    for (int i = 0; i < size; ++i) {
        temperatures[i] = start + i * step;
    }

    // Stampiamo le temperatures per verificarne il contenuto
    for (int i = 0; i < size; ++i) {
        std::cout << temperatures[i] << " ";
    }


  
  System SYS;
  
  cout<<"parte simulazione" << endl;
  
  //////SIMULATION////////
  for(int l=0; l < size ; l++){
      SYS.initialize_for_iteration(temperatures[l]);
      //SYS.set_temperature(temperatures[l]);
      SYS.initialize_properties();
      if (l==0) {SYS.open_file_T();}
      SYS.block_reset(0);
      cout << SYS.get_temp() << "   " << temperatures[l] << endl;
      for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
          for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
            SYS.step();
            SYS.measure();
          }
          cout<< "medie blocco "<< i+1 << endl;
          SYS.averages(i+1);
          SYS.block_reset(i+1);
        }
      SYS.finalize();}




  cout<<"programma eseguito"<<endl;
  return 0;
}

