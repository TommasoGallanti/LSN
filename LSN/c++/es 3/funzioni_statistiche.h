#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>
//#include "TH1F.h"
//#include "TApplication.h"
//#include "TCanvas.h"

using namespace std ;

 template <typename> vector<double> LeggiDatiFile (unsigned int N , string nome_file ) {

  vector <double> v;
  double value = 0  ;
  ifstream fin (nome_file);

  if (!fin) {
       cout << " empty " << endl ; 
       exit(0);
    }
  else {
      for (unsigned int k = 0 ; k < N ; k++ ) {
          
          fin >> value ;
          v.push_back(value) ;

          if ( fin.eof()) {
              cout << " end of file " << endl;
              exit (0) ;
            }     
        }
    }
 return v ;
 }  



 template <typename T> double CalcolaMedia ( const vector <T> & v ) {

 double media = 0 ; 

 if (v.size()==0) return media;

 
 for (int k=0; k < v.size() ; k++ ) {
     media = media + v[k] ; 
     } 
      
 return media/double(v.size()); 
 }  


 template <typename> double CalcolaVarianza ( const vector<double> & v , int media ) {
    
    double varianza = 0 ; 
 
    if (v.size()==0) return varianza;

     for (int k=0; k< v.size() ; k++ ) {
        varianza = varianza + ((media - v[k]) * (media - v[k])) ;
        }
     
     return varianza/double(v.size()) ;
     }



 template <typename> vector<double> Copiavector (vector <double> & v ) {
    vector<double> n (v.size());
    
    for (int i = 0 ; i < v.size() ; i++)  n[i] = v[i];
    return n ;
}

 template <typename> void swap ( unsigned int a , unsigned int b, vector<double> &  v) {
    
    double t= v[a] ;
    v[a] = v[b] ; 
    v[b] = t ;
    
}




 template <typename> void SelSort (vector <double> & v ) {
        for ( int k = 0 ; k< v.size() ; k++) { 
             for ( int j = k+1 ; j < v.size(); j++ ) { 
                 if (v[j] > v[k]) { 
                    swap<void> (j,k,v) ; 
                }   
            }
        }    
    
}


 template <typename> double CalcolaMediana (vector <double> &v) {
    
    SelSort<void> (v) ;
    double mediana = 0; 
    
 if (v.size()%2==0) {  mediana = (v[v.size()/2] + v[v.size()/2-1])/2 ; }
 
 else { mediana = v[double(v.size())/2]  ;  }
 
 return mediana;
}


template <typename> void print (vector<double> &v , string filename) {
    
    SelSort<void> (v) ;
    ofstream fout (filename) ;
  for ( int k = 0 ; k < v.size() ; k++ ) fout << v[k] << endl; 
  fout.close();
  return;
}