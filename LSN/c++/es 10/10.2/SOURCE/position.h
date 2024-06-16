#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "random.h"
#include <vector>

using namespace std;
using namespace arma;

#pragma once


class Position {
private:
    double _x;
    double _y;
    
    
    
    // Random number generator
    Random _rnd;

public:
    // Constructor
    //Position() : _x(0), _y(0), _rnd() {}
    
    // Getter methods
    double get_x();
    double get_y();
    double get_r();

    // Setter methods
    void set_x(double x);
    void set_y(double y);
    void Initialize();

    // Utility methods
    static double get_distance(Position p1, Position p2);
    void print_xy();
    void print_xy_file();
    std::vector<Position> create_map(int N_city, int _sim_type, double r = 1.0);
    std::vector<Position> read_map(string filename);
    int countLines(const string& filename);

    // Optional: method to check if the position fulfills certain bounds
    // bool check_bounds(double x_min, double x_max, double y_min, double y_max);
};


