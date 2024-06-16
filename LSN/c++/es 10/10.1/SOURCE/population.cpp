#include <cmath>
#include <cstdlib>
#include <string>
#include "population.h"
#include "position.h"
#include <armadillo>
#include <vector>
#include <set>
#include <unordered_map>


using namespace std;
using namespace arma;






void Population::Initialize(int rank){
    int p1, p2; // Read from Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);
  ofstream coutf;
  coutf.open("initial_population.dat");
  coutf.close();
  coutf.open("check_population.dat");
  coutf.close();
  coutf.open("Mean_loss_"+to_string(rank)+".dat");
  coutf << "# generation:     mean_loss:           best_loss:" << endl;
  coutf.close();
  Position pos;
  int n_city = pos.countLines("../INPUT/cap_prov_ita.dat");
  this->set_n_city(n_city);
}


vector<int> Population::riordina_array(const vector<int>& primo_array, const vector<int>& secondo_array) {
    // Crea una mappa che associa ciascun numero nel secondo array con il suo indice
    unordered_map<int, int> indice_numeri;
    for (int i = 0; i < secondo_array.size(); ++i) {
        indice_numeri[secondo_array[i]] = i;
    }
    
    // Ordina il primo array basandosi sull'ordine dei numeri nel secondo array
    vector<int> risultato(primo_array);
    sort(risultato.begin(), risultato.end(), [&](int a, int b) {
        return indice_numeri[a] < indice_numeri[b];
    });
    
    return risultato;
}


void Population::sortAscending(std::vector<double>& vec) {
    // Utilizziamo std::sort per ordinare in ordine crescente
    std::sort(vec.begin(), vec.end());
}


void Population::shuffle(std::vector<int>& array, int startIndex, int endIndex) {
    if (endIndex >= array.size()-1) {endIndex = array.size() -2;}
    std::reverse(array.begin() + startIndex, array.begin() + endIndex + 1);
}

void Population::single_mutation2(std::vector<int>& array, int startIndex) {
    int m = int (_rnd.Rannyu(1.,_n_city-1));
    int n = 2;
    int endIndex = startIndex + m;
    if (endIndex >= array.size() - 1) {
        endIndex = array.size() - 2; // Imposta l'indice finale all'ultimo elemento meno uno
    }
    std::vector<int> subarray(array.begin() + startIndex, array.begin() + endIndex + 1);
    std::rotate(subarray.rbegin(), subarray.rbegin() + n, subarray.rend());
    std::copy(subarray.begin(), subarray.end(), array.begin() + startIndex);
}

void Population::single_mutation3(std::vector<int>& array, int startIndex1, int startIndex2) {
    int m = 3;
    int endIndex1 = startIndex1 + m ;
    int endIndex2 = startIndex2 + m ;
    if (endIndex1 >= array.size() - 1 || endIndex2 >= array.size() - 1) {
        endIndex1 = std::min(endIndex1, static_cast<int>(array.size()) - 2); // Imposta l'indice finale all'ultimo elemento meno uno
        endIndex2 = std::min(endIndex2, static_cast<int>(array.size()) - 2); // Imposta l'indice finale all'ultimo elemento meno uno
    }

    // Permuta le due sottosequenze
    for (int i = 0; i < m; ++i) {
        std::swap(array[startIndex1 + i], array[startIndex2 + i]);
    }
}

std::pair<std::vector<int>, std::vector<int>> Population::Crossover(const std::vector<int>& mother, const std::vector<int>& father) {
    int n = mother.size();
    // Seleziona casualmente il punto di taglio all'interno del percorso, escludendo il primo e l'ultimo elemento
    int crossoverPoint = int((_n_city-3)*_rnd.Rannyu())+1;
    //cout << crossoverPoint << endl;
    std::vector<int> child1;
    std::vector<int> child2;

    // Conserva le prime parti dei percorsi dei genitori nei figli
    for (int i = 0; i < crossoverPoint; ++i) {
        child1.push_back(mother[i]);
        child2.push_back(father[i]);
    }

    

    // Completa i percorsi dei figli con le città mancanti dall'altro genitore
    int indexChild1 = crossoverPoint;
    int indexChild2 = crossoverPoint;
    std::vector<int> city_missing1;
    std::vector<int> city_missing2;
    std::vector<int> city_missing1_ordinate;
    std::vector<int> city_missing2_ordinate;
    for (int i = crossoverPoint; i < n-1; ++i) {
        //inserisci in city missing1 le città mancanti 
        city_missing1.push_back(mother[i]); 

        //inserisci in city missing 2 le città mancanti 
        city_missing2.push_back(father[i]);
    }
    city_missing1_ordinate = riordina_array(city_missing1,father); //ordina le città mancanti
    city_missing2_ordinate = riordina_array(city_missing2,mother);

    for (int i = 0; i < city_missing1_ordinate.size(); ++i) {
        //inserisci le città mancanti ordinate nei child
        child1.push_back(city_missing1_ordinate[i]);
        child2.push_back(city_missing2_ordinate[i]);
    }

    child1.push_back(1);
    child2.push_back(1);

    return std::make_pair(child1, child2);
}




void Population::sortIndividualsByLoss(vector<Individuo_Loss>& data) {
    std::sort(data.begin(), data.end(), [](const Individuo_Loss& a, const Individuo_Loss& b) {
        return a.loss < b.loss; // Ordine crescente di loss
    });
}

void Population::printIndividualsWithLosses(const vector<Individuo_Loss>& data) {
    for (const auto& item : data) {
        cout << "Indice originale: " << item.original_index << " - Loss: " << item.loss << endl;
    }
}

void Population::mutation(std::vector<int>&  individuo) {
    double pm1 = 0.25;
    double pm2 = 0.5;
    double pm3 = 0.75;
    double pm4 = 1.;
    double extraction = _rnd.Rannyu();
    //cout << extraction << endl;
    int index1 = int((_n_city-3)*_rnd.Rannyu())+1;
    int index2 = int((_n_city-3)*_rnd.Rannyu())+1;
    if (extraction<=pm1) {
        std::swap(individuo[index1],individuo[index2]);
        mut1++;
    }            //scambia due città
    else if (pm1 < extraction && extraction <= pm2) {
        mut2++;
        single_mutation2(individuo,index1);
        //cout << "mut2" << endl;
    }        //shift of +n positions
    else if (pm2 < extraction && extraction <= pm3) {
        single_mutation3(individuo,index1,index2);
        mut3++;
    }        //scambia seconda con penultima e terza con terzultima
    else if (pm3 < extraction && extraction <= pm4) {
        mut4++;
        if (index1<index2) {
            shuffle(individuo,index1,index2);
        }
        else if (index1>index2) {
            shuffle(individuo,index2,index1);
        }
        else {}
    }        //shuffle
    
}

std::vector<std::vector<int>> Population::create_evolution(std::vector<Individuo_Loss> individuals_with_losses){
    // ATTENZIONE!!! IL VETTORE CONTENENTE GLI STRUCT DEVE ESSERE ORDINATO IN FUNZIONE DELLA LOSS
    // 0 - CREATE VECTOR FOR EVO_POPULATION
    ofstream out;
    out.open("indici.dat");
    double p_cross = 0.8;
    double p_mutation = 0.2;
    std::vector<std::vector<int>> evo_population; // Popolazione evoluta
    evo_population.reserve(individuals_with_losses.size()); // Riserva memoria per migliorare le prestazioni
    
    while (evo_population.size() < individuals_with_losses.size()) {   // Ripeti fino a quando la popolazione evoluta ha la stessa dimensione della popolazione iniziale
        
        // 1 - SELECT A PAIR OF INDIVIDUALS IN FUNCTION OF LOSS
        int index1 = int(_n_individui*std::pow(_rnd.Rannyu(),5.5));
        int index2 = int(_n_individui*std::pow(_rnd.Rannyu(),5.5));
        out << "indici: " << index1 << "  "<< index2 <<  endl;
        std::vector<int> ind1 = individuals_with_losses[index1].individual;
        std::vector<int> ind2 = individuals_with_losses[index2].individual;
        //std::vector<int> child1 = individuals_with_losses[index1].individual;
        //std::vector<int> child2 = individuals_with_losses[index2].individual;
        std::vector<int> child1;
        std::vector<int> child2;
        // 2 - DO A CROSSOVER BETWEEN THE TWO INDIVIDUALS SELCTED WITH PROBABILITY P_crossover
        if (_rnd.Rannyu()<p_cross) {
            std::pair<std::vector<int>, std::vector<int>> children = Crossover(ind1, ind2);
            child1 = children.first;
            child2 = children.second;
        }
        else {
            child1 = individuals_with_losses[index1].individual;
            child2 = individuals_with_losses[index2].individual;
        }
        // 3 - DO A MUTATION OF THE INDIVIDUALS WITH PROBABILITY P_mutation
        if (_rnd.Rannyu()<p_mutation) {
            mutation(child1);
            mutation(child2);
        }
        else {}

        // 4 - ADD THE TWO NEW INDIVIDUALS AT THE EVO_POPULATION
        evo_population.push_back(child1);
        if (evo_population.size() < individuals_with_losses.size()) {
            evo_population.push_back(child2);
        }
    }   // 5 - REPEAT FROM STEP 2 TO STEP 5 UNTIL THE EVO_POPULATION HAS THE SAME DIMENSION OF THE INITIAL POPULATION 
    out.close();
    return evo_population; 
}


std::vector<std::vector<int>> Population::create_population(const int N_city, int n_individui) {
    std::vector<std::vector<int>> population; // Popolazione di vettori di interi
    std::vector<int> new_pop(N_city-1); // Dimensione N_city (senza l'elemento extra per il valore iniziale e finale)
    std::iota(new_pop.begin(), new_pop.end(), 2); // Riempie con valori da 2 a N_city

    population.reserve(n_individui); // Riserva memoria per migliorare le prestazioni

    // Inizializza il generatore di numeri casuali
    std::random_device rd;
    std::mt19937 g(rd());

    // Genera n_individui popolazioni casuali
    for (int i = 0; i < n_individui; ++i) {
        // Mescola l'ordine dei numeri in new_pop
        std::shuffle(new_pop.begin(), new_pop.end(), g);

        // Crea un nuovo vettore che inizia e finisce con 1
        std::vector<int> individual(N_city + 1);
        individual[0] = 1;
        std::copy(new_pop.begin(), new_pop.end(), individual.begin() + 1);
        individual[N_city] = 1;

        // Inserisce il nuovo vettore mescolato nella popolazione
        population.push_back(individual);
    }

    return population;
}

void Population::print_on_file(std::vector<std::vector<int>> population){
    ofstream coutf;
    std::vector<int> individuo;
    coutf.open("initial_population.dat",ios::app);
    for (int i=0; i<_n_individui; i++){
        //individuo = this->get_individuo(i,population);
        individuo = population[i];
        for (int l=0; l<_n_city+1; l++){
            //individuo = this->get_individuo(i,population);
            if (l==_n_city){coutf << individuo[l] << endl;}
            else {coutf << individuo[l] << setw(15);} 
        }
    }
}

std::vector<int> Population::get_individuo(int index_individuo, std::vector<std::vector<int>> population){
    return population[index_individuo];
}

bool Population::check_population(std::vector<std::vector<int>> population){
    ofstream couta;
    couta.open("check_population.dat", ios::app);
    bool _first_value = true;
    bool _last_value = true;
    bool no_rep=true;
    bool check = false;
    std::set<int> unique_elements; // Set per tenere traccia degli elementi univoci
    for (int i=0; i<_n_individui; i++){
        //std::vector<int> individuo = this->get_individuo(i,population);
        std::vector<int> individuo = population[i];
        if (individuo[0]!=1) {_first_value=false;}
        if (individuo[_n_city]!=1) {_last_value=false;}
        unique_elements.clear(); // Resetta il set per ogni individuo
        for (int l=1; l<_n_city; l++){
            //voglio che tutti gli elementi dentro individuo siano ripetuti solo una volta
            if (unique_elements.find(individuo[l]) != unique_elements.end()) {
                // Se l'elemento è già presente, setta no_rep a false e interrompi il loop
                no_rep = false;
                break;
            } else {
                // Altrimenti, aggiungi l'elemento al set
                unique_elements.insert(individuo[l]);
            }
        }
        if (!no_rep) {
            // Se uno degli individui ha elementi duplicati, esci dal loop
            break;
        }
        
        //couta << "individuo " << i+1 << " = " << _first_value << "  " << _last_value << "  " << no_rep << endl;
    }
    couta.close();
    if (_first_value == true && _last_value == true && no_rep && true) {check = true;}
    return check;
}

double Population::get_loss(std::vector<Position> trip){
    double loss=0 ;
    for (int i =0; i < _n_city; i++){
        loss += Position::get_distance(trip[i],trip[i+1]);  
    }
    //loss += Position::get_distance(trip[_n_city],trip[0]);
    return loss;
}

std::vector<Position> Population::create_trip(std::vector<Position> map, std::vector<int> individuo){
    std::vector<Position> trip;
    trip.reserve(_n_city+1); // Riserva memoria per migliorare le prestazioni
    int index;
    for(int i=0; i < _n_city+1; i++){
        index = individuo[i]-1;
        trip.push_back(map[index]);
    }
    return trip;
}

void Population::print_loss(std::vector<double> loss){
    ofstream coutf;
    coutf.open("Loss.dat");
    for (int i=0; i<loss.size(); i++){
        //coutf << loss[i] << endl;
    }
}

void Population::print_mutation_rate(){
    ofstream out;
    out.open("mutations");
    out << mut1 << "  " << mut2 << "  " << mut3 << "  "<< mut4 << endl;
}



void Population::set_n_city(int n_city) {
    _n_city = n_city;
}
