#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include "random.h"
#include "position.h"

class Population {
private:
    Random _rnd;  // Random number generator
    int _n_individui = 1000;
    int _n_city = 34;
    void shuffle(std::vector<int>& array, int startIndex, int endIndex);
    void single_mutation2(std::vector<int>& array, int startIndex);
    void single_mutation3(std::vector<int>& array, int startIndex1, int startIndex2);
    void mutation(std::vector<int>& individuo);
    std::pair<std::vector<int>, std::vector<int>> Crossover(const std::vector<int>& mother, const std::vector<int>& father);
    int mut1;
    int mut2;
    int mut3;
    int mut4;

public:
    

    // Initializes the random number generator
    void Initialize();

    // Creates a population of size n_individui with cities numbered from 1 to N_city
    std::vector<std::vector<int>> create_population(const int N_city, int n_individui);

    // Prints the population to a file
    void print_on_file(std::vector<std::vector<int>> population);
    void print_loss(std::vector<double> loss);
    bool check_population(std::vector<std::vector<int>> population);
    void print_mutation_rate();

    // Returns an individual (vector of integers) from the population at a given index
    std::vector<int> get_individuo(int index_individuo, std::vector<std::vector<int>> population);

    double get_loss(std::vector<Position> trip);

    std::vector<Position> create_trip(std::vector<Position> map, std::vector<int> individuo);
    
    //generation
    struct Individuo_Loss {
        std::vector<int> individual;
        double loss;
        int original_index;

        // Costruttore per inizializzare la struttura
        Individuo_Loss(std::vector<int> ind, double l, int idx) : individual(ind), loss(l), original_index(idx) {}
    };

    void sortAscending(std::vector<double>& vec);
    void sortIndividualsByLoss(std::vector<Individuo_Loss>& data);
    void printIndividualsWithLosses(const std::vector<Individuo_Loss>& data);
    std::vector<std::vector<int>> create_evolution(std::vector<Individuo_Loss> individuals_with_losses);
    vector<int> riordina_array(const vector<int>& primo_array, const vector<int>& secondo_array);

};


