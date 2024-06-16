#include <iostream>
#include "position.h"
#include "population.h"
#include <armadillo>
#include "random.h"
#include <vector>
#include <algorithm> // Necessario per std::sort


using namespace std;
using namespace arma;

Random _rnd;

void sortAscending(std::vector<double>& vec) {
    // Utilizziamo std::sort per ordinare in ordine crescente
    std::sort(vec.begin(), vec.end());
}

struct Individuo_Loss {
    vector<int> individual;
    double loss;
    int original_index;  // Aggiungi questo campo per tracciare l'indice originale

    // Costruttore per inizializzare la struttura
    Individuo_Loss(vector<int> ind, double l, int idx) : individual(ind), loss(l), original_index(idx) {}
};

void sortIndividualsByLoss(vector<Individuo_Loss>& data) {
    std::sort(data.begin(), data.end(), [](const Individuo_Loss& a, const Individuo_Loss& b) {
        return a.loss < b.loss; // Ordine crescente di loss
    });
}

void printIndividualsWithLosses(const vector<Individuo_Loss>& data) {
    for (const auto& item : data) {
        cout << "Indice originale: " << item.original_index << " - Loss: " << item.loss << endl;
    }
}



void shuffle(std::vector<int>& array, int startIndex, int endIndex) {
    std::reverse(array.begin() + startIndex, array.begin() + endIndex + 1);
}

void single_mutation2(std::vector<int>& array, int startIndex, int m, int n) {
    int endIndex = startIndex + m;
    std::vector<int> subarray(array.begin() + startIndex, array.begin() + endIndex + 1);
    std::rotate(subarray.rbegin(), subarray.rbegin() + n, subarray.rend());
    std::copy(subarray.begin(), subarray.end(), array.begin() + startIndex);
}

void single_mutation3(std::vector<int>& array, int startIndex1, int m, int startIndex2) {
    int endIndex1 = startIndex1 + m ;
    int endIndex2 = startIndex2 + m ;

    // Permuta le due sottosequenze
    for (int i = 0; i < m; ++i) {
        std::swap(array[startIndex1 + i], array[startIndex2 + i]);
    }
}

void mutation(std::vector<int> individuo,const int _n_city) {
    double pm1 = 0.09;
    double pm2 = 0.18;
    double pm3 = 0.27;
    double pm4 = 0.36;
    double extraction = _rnd.Rannyu();
    int index1 = int((_n_city-2)*_rnd.Rannyu())+1;
    int index2 = int((_n_city-2)*_rnd.Rannyu())+1;
    if (extraction<=pm1) {std::swap(individuo[index1],individuo[index2]);}            //scambia due città
    else if (pm1<extraction<=pm2) {single_mutation2(individuo,index1,4,2);}        //shift of +n positions
    else if (pm2<extraction<=pm3) {single_mutation3(individuo,index1,6,index2);}        //scambia seconda con penultima e terza con terzultima
    else if (pm3<extraction<=pm4) {
        if (index1<index2) {
            shuffle(individuo,index1,index2);
        }
        else if (index1>index2) {
            shuffle(individuo,index2,index1);
        }
    }        //shuffle
}

std::vector<std::vector<int>> create_evolution(std::vector<Individuo_Loss> individuals_with_losses,const int _n_individui, const int _n_city){
    // ATTENZIONE!!! IL VETTORE CONTENENTE GLI STRUCT DEVE ESSERE ORDINATO IN FUNZIONE DELLA LOSS
    // 0 - CREATE VECTOR FOR EVO_POPULATION
    double p_cross = 0.6;
    std::vector<std::vector<int>> evo_population; // Popolazione evoluta
    evo_population.reserve(individuals_with_losses.size()); // Riserva memoria per migliorare le prestazioni
    
    while (evo_population.size() < individuals_with_losses.size()) {   // Ripeti fino a quando la popolazione evoluta ha la stessa dimensione della popolazione iniziale
        
        // 1 - SELECT A PAIR OF INDIVIDUALS IN FUNCTION OF LOSS
        int index1 = int(_n_individui*pow(_rnd.Rannyu(),10));
        int index2 = int(_n_individui*pow(_rnd.Rannyu(),10));
        std::vector<int> ind1 = individuals_with_losses[index1].individual;
        std::vector<int> ind2 = individuals_with_losses[index2].individual;

        // 2 - DO A CROSSOVER BETWEEN THE TWO INDIVIDUALS SELCTED WITH PROBABILITY P_crossover
        if (_rnd.Rannyu() < p_cross) {
            int index_cross1 = int((_n_city-2)*_rnd.Rannyu())+1;
            int index_cross2 = int((_n_city-2)*_rnd.Rannyu())+1;
            if (index_cross1<index_cross2) {
                shuffle(ind1,index_cross1,index_cross2);
                shuffle(ind2,index_cross1,index_cross2);
            }
            else if (index_cross1>index_cross2) {
                shuffle(ind1,index_cross2,index_cross1);
                shuffle(ind2,index_cross2,index_cross1);
            }
        }

        // 3 - DO A MUTATION OF THE INDIVIDUALS WITH PROBABILITY P_mutation
        mutation(ind1,_n_city);
        mutation(ind2,_n_city);

        // 4 - ADD THE TWO NEW INDIVIDUALS AT THE EVO_POPULATION
        evo_population.push_back(ind1);
        if (evo_population.size() < individuals_with_losses.size()) {
            evo_population.push_back(ind2);
        }
    }   // 5 - REPEAT FROM STEP 2 TO STEP 5 UNTIL THE EVO_POPULATION HAS THE SAME DIMENSION OF THE INITIAL POPULATION    
    return evo_population; 
}




int main (int argc, char *argv[]){
    
    int n_individui = 1000;
    int N_city = 34;
    int _sim_type=0;
    cout << "numero città = " << N_city << endl;
    cout << "numero individui = " << n_individui << endl;
    if (_sim_type==0) {cout << "mappa circolare" << endl;}
    else if (_sim_type==1) {cout << "mappa quadrata" << endl;}
    else {
        cout << "error! invalid sym_type, ammessi solo 0 o 1";
        exit(EXIT_FAILURE);
    }
    
    double r=1.;
    Position pos;
    pos.Initialize();

    vector <Position> positions=pos.create_map(N_city,_sim_type,r);
    cout << "mappa creata" << endl;

    Population population;
    population.Initialize();

    vector <vector<int>> population_0 = population.create_population(N_city,n_individui);
    population.print_on_file(population_0);
    cout << "popolazione iniziale creata" << endl;
    
    bool check = population.check_population(population_0);
    if (check==true) {cout << "check su population passato!" << endl;}
    else if (check==false) {cout << "check su population non passato!" << endl;}
    
    //cout << population_0[0][0] << endl;
    //cout << Position::get_distance(positions[0],positions[1]) << endl;

    vector<double> losses;
    vector<Individuo_Loss> individui_with_losses;
    for (int i=0; i<n_individui; i++){
        std::vector<Position> trip = population.create_trip(positions,population_0[i]);
        losses.push_back(population.get_loss(trip));
        individui_with_losses.emplace_back(population_0[i],population.get_loss(trip),i); //creo lo struct con, individuo, loss e indice individuo
    }
    
    //sortAscending(losses);
    sortIndividualsByLoss(individui_with_losses); //ordino gli struct in base alla loss 
    cout << "Individui con Loss ordinati:" << endl;
    //printIndividualsWithLosses(individui_with_losses); //stamnpo su terminale loss e index_individuo

    // Stampo info su best individuo
    const Individuo_Loss& best_individual = individui_with_losses[0];
    cout << "Individuo con la loss più bassa:" << endl;  
    cout << "Indice originale: " << best_individual.original_index << endl;
    cout << "Loss: " << best_individual.loss << endl;
    cout << "Individuo per esteso: ";
    for (int gene : best_individual.individual) {
        cout << gene << " ";
    }
    cout << endl;

    vector <vector<int>> population_1 = create_evolution(individui_with_losses,n_individui,N_city);
    population.print_on_file(population_1);
    cout << "prima generazione creata" << endl;
    
    check = population.check_population(population_1);
    if (check==true) {cout << "check su population1 passato!" << endl;}
    else if (check==false) {cout << "check su population1 non passato!" << endl;}

    ///////
    vector<double> losses1;
    vector<Individuo_Loss> individui_with_losses1;
    for (int i=0; i<n_individui; i++){
        std::vector<Position> trip = population.create_trip(positions,population_1[i]);
        losses1.push_back(population.get_loss(trip));
        individui_with_losses1.emplace_back(population_1[i],population.get_loss(trip),i); //creo lo struct con, individuo, loss e indice individuo
    }
    
    //sortAscending(losses);
    sortIndividualsByLoss(individui_with_losses1); //ordino gli struct in base alla loss 
    cout << "Individui con Loss ordinati:" << endl;
    //printIndividualsWithLosses(individui_with_losses); //stamnpo su terminale loss e index_individuo

    // Stampo info su best individuo
    const Individuo_Loss& best_individual1 = individui_with_losses1[0];
    cout << "Individuo con la loss più bassa:" << endl;  
    cout << "Indice originale: " << best_individual1.original_index << endl;
    cout << "Loss: " << best_individual1.loss << endl;
    cout << "Individuo per esteso: ";
    for (int gene : best_individual1.individual) {
        cout << gene << " ";
    }
    cout << endl;

    return 0;


}
