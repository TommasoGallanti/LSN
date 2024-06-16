#include <cmath>
#include <cstdlib>
#include <string>
#include "position.h"
#include <armadillo>

using namespace std;
using namespace arma;



void Position::Initialize(){
    int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);
  ofstream coutf;
  coutf.open("map.dat");
  coutf.close();

}


double Position::get_x(){
    return _x;
}

double Position::get_y(){
    return _y;
}

double Position::get_r(){
    return sqrt(pow(_x,2)+pow(_y,2));
}

void Position::set_x(double x){
    _x=x;
}

void Position::set_y(double y){
    _y=y;
}

double Position::get_distance(Position p1, Position p2){
    return pow(p1.get_x() - p2.get_x(),2) + pow(p1.get_y() - p2.get_y(),2);
}

void Position::print_xy(){
    cout << "x = " << this->get_x() << "  y = " << this->get_y() << endl ;
}

void Position::print_xy_file(){
    ofstream couta;
    couta.open("map.dat",ios::app);
    couta << this->get_x() << "         " << this->get_y() << endl ;
    couta.close();
}



//bool Position::check_bonds(){
  //  if (){return true}
  //  else {return false}
//}

std::vector<Position> Position::create_map(int N_city, int _sim_type, double r) {
    //Random rnd;
    //rnd.StartGen();
    std::vector<Position> positions;
    positions.reserve(N_city); // Riserva memoria per migliorare le prestazioni
    double theta = 0.;
    for (int i = 0; i < N_city; ++i) {
        Position new_position; // Crea un nuovo oggetto Position
        //new_position.Initialize();
        if (_sim_type == 1) {
            new_position.set_x(_rnd.Rannyu());
            new_position.set_y(_rnd.Rannyu());
            //new_position.print_xy();
            new_position.print_xy_file();
            
        } else if (_sim_type == 0) {
            theta = _rnd.RanAngle();
            //cout << "iterazione " << i << " theta = "<< theta << endl;
            new_position.set_x(r * cos(theta));
            new_position.set_y(r * sin(theta));
            //new_position.print_xy();
            new_position.print_xy_file();
        }

        positions.push_back(new_position); // Inserisce il nuovo oggetto nel vettore
    }

    return positions;
}


