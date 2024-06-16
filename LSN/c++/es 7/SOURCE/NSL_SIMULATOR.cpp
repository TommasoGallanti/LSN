/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include "system.h"
#include <armadillo>

using namespace std;
using namespace arma;

int main (int argc, char *argv[]){

  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);
  
  if (SYS.get_simtype()==0){    ////effettua equilibration in caso di NVE
    cout<< "Fase di equilibrazione partendo da T = " << SYS.get_temp() << endl ;
    cout << "..." << endl;
    ///////EQUILIBRATION//////
    for(int i=0; i < SYS.get_nbl_eq(); i++){ //loop over blocks
        for(int j=0; j < SYS.get_nsteps_eq(); j++){ //loop over steps in a block
          SYS.step();
          SYS.measure();
        }
      SYS.averages_eq(i+1);
      SYS.block_reset(i+1);
    }
    
    cout<<"dopo " << SYS.get_nbl_eq() << " blocchi da " << SYS.get_nsteps_eq() << " steps, viene raggiunto l'equilibrio, per l'ultimo step" << endl;
    cout<<"T_eq = " << SYS.get_temp_equilibration() << endl;
    SYS.block_reset(0);
    SYS.set_temp();    //setto temperatura dopo eq
    //cout<<"setto temperatura a " << SYS.get_temp() << endl;

    cout<< "Una volta raggiunto l'equilibrio, setto la temperatura desiderata a " << SYS.get_temp() << endl << "e riavvio la simulazione" << endl;
    cout<< "numero blocchi " << SYS.get_nbl() << endl;
    cout<< "numero step " << SYS.get_nsteps() << endl;
  }
  
  
  //////SIMULATION////////
  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
      }
      cout<< "medie blocco "<< i+1 << endl;
      SYS.averages(i+1);
      SYS.block_reset(i+1);
    }
  SYS.finalize();
  


  int n_attemps = SYS.get_nattempts();
  int n_accepted = SYS.get_naccepted();
  cout << "n_attemps = "<< SYS.get_nattempts() << endl;
  cout << "n_accepted = " << SYS.get_naccepted() << endl;
  cout << "acceptance ratio = " << double(n_accepted)/double(n_attemps) << endl;

  //if (SYS.get_gofr_permission()==1) SYS.print_gofr();    /////ATTENZIONE!!! 
  cout<<"programma eseguito"<<endl;
  return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
