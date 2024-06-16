#include <iostream>
#include "position.h"
#include "population.h"
#include <armadillo>
#include "random.h"
#include <vector>
#include <algorithm> // Necessario per std::sort


using namespace std;
using namespace arma;


double calcolaMedia_first_half(const std::vector<double>& vettore) {
    double somma = 0.0;
    int dimensione_metà = vettore.size() / 2; // Dimensione della prima metà del vettore

    // Somma solo gli elementi della prima metà del vettore
    for (int i = 0; i < dimensione_metà; ++i) {
        somma += vettore[i];
    }
    
    // Calcola la media della prima metà del vettore
    return somma / dimensione_metà;
}
int main (int argc, char *argv[]){
    
    const int n_individui = 1000;
    const int N_city = 34;
    int _sim_type=1;
    int n_generations = 200;
    cout << "numero città = " << N_city << endl;
    cout << "numero individui = " << n_individui << endl;
    cout << "numero generazioni = " << n_generations << endl; 
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
    if (check==true) {cout << "check su population iniziale passato!" << endl;}
    else if (check==false) {cout << "check su population non passato!" << endl;}
    
    //cout << population_0[0][0] << endl;
    //cout << Position::get_distance(positions[0],positions[1]) << endl;

    vector<double> losses;
    vector<Population::Individuo_Loss> individui_with_losses;
    for (int i=0; i<n_individui; i++){
        std::vector<Position> trip = population.create_trip(positions,population_0[i]);
        if (i==0) {cout << trip.size() << endl;}
        losses.push_back(population.get_loss(trip));
        individui_with_losses.emplace_back(population_0[i],population.get_loss(trip),i); //creo lo struct con, individuo, loss e indice individuo
    }
    
    //sortAscending(losses);
    population.sortIndividualsByLoss(individui_with_losses); //ordino gli struct in base alla loss 
    //cout << "Individui con Loss ordinati:" << endl;
    //printIndividualsWithLosses(individui_with_losses); //stamnpo su terminale loss e index_individuo
    ofstream out;
    out.open("Mean_loss.dat",ios::app);
    for (int generation = 0; generation < n_generations; generation++) {
        vector<vector<int>> population_1 = population.create_evolution(individui_with_losses);
        population.print_on_file(population_1);
        

        check = population.check_population(population_1);
        if (check == true) {} //cout << "check passato!" << endl; }
        else if (check == false) { 
            cout << "check su population " << generation + 1 << " non passato!" << endl; 
            break; 
        }

        vector<double> losses1;
        vector<Population::Individuo_Loss> individui_with_losses1;
        for (int i = 0; i < n_individui; i++) {
            vector<Position> trip = population.create_trip(positions, population_1[i]);
            double loss = population.get_loss(trip);
            losses1.push_back(loss);
            individui_with_losses1.emplace_back(population_1[i], loss, i);
        }

        population.sortIndividualsByLoss(individui_with_losses1);
        population.sortAscending(losses1);
        population.print_loss(losses1);
        out << generation+1 << setw(25) << calcolaMedia_first_half(losses1) << setw(25) << individui_with_losses1[0].loss << endl;
        
        //cout << "Generazione " << generation + 1 << endl;
        //const Population::Individuo_Loss best_individual1 = individui_with_losses1[0];
        //cout << "Individuo con la loss più bassa:" << endl;  
        //cout << "Indice originale: " << individui_with_losses1[0].original_index << endl;
        cout << "best Loss: " << individui_with_losses1[0].loss << endl;
        

        // Aggiorna individui_with_losses per la prossima generazione
        individui_with_losses = individui_with_losses1;
    }
    ofstream outf;
    outf.open("final_trip.dat");
    vector<Position> final_trip = population.create_trip(positions, individui_with_losses[0].individual);
    for (int i =0; i<final_trip.size();i++){
        outf << final_trip[i].get_x() << "     " << final_trip[i].get_y() << endl;
    }
    outf.close();
    population.print_mutation_rate();
    return 0;


}
