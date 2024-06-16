#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "position.h"
#include "population.h"
#include <armadillo>
#include "random.h"
#include <vector>
#include <algorithm>
#include <mpi.h>

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

int main(int argc, char* argv[]) {

    
    int size ; 
    int rank ;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int numero_change =0;
    
    const int n_individui = 1000;
    int n_generations = 2000;
    cout << "numero individui = " << n_individui << endl;
    
    Position pos;
    pos.Initialize();

    vector <Position> positions = pos.read_map("../INPUT/cap_prov_ita.dat");
    cout << "mappa letta" << endl;
    int N_city = pos.countLines("../INPUT/cap_prov_ita.dat");
    cout << "numero città = " << N_city << endl;

    Population population;
    population.Initialize(rank + 1);
    

    vector <vector<int>> population_0 = population.create_population(N_city, n_individui);
    population.print_on_file(population_0);
    cout << "popolazione iniziale creata" << endl;

    std::vector<Position> trip = population.create_trip(positions, population_0[0]);
    double best_loss = population.get_loss(trip);
    
    bool check = population.check_population(population_0);
    if (check == true) {
        cout << "check su population iniziale passato!" << endl;
    } else if (check == false) {
        cout << "check su population non passato!" << endl;
    }

    vector<double> losses;
    vector<Population::Individuo_Loss> individui_with_losses;
    for (int i = 0; i < n_individui; i++) {
        std::vector<Position> trip = population.create_trip(positions, population_0[i]);
        if (i == 0) {
            cout << "trip size = " << trip.size() << endl;
        }
        losses.push_back(population.get_loss(trip));
        individui_with_losses.emplace_back(population_0[i], population.get_loss(trip), i); //creo lo struct con, individuo, loss e indice individuo
    }
    
    population.sortIndividualsByLoss(individui_with_losses); //ordino gli struct in base alla loss 

    ofstream out;
    out.open("Mean_loss_" + to_string(rank + 1) + ".dat", ios::app);

    vector<double> x_trans(N_city);
    vector<double> y_trans(N_city);

    for (int generation = 0; generation < n_generations; generation++) {
        vector<vector<int>> population_1 = population.create_evolution(individui_with_losses);
        
        check = population.check_population(population_1);
        if (check == false) {
            cout << "check su population " << generation + 1 << " non passato!" << endl;
            break;
        }

        vector<Position> trip;
        vector<double> losses1;
        vector<Population::Individuo_Loss> individui_with_losses1;
        for (int i = 0; i < n_individui; i++) {
            trip = population.create_trip(positions, population_1[i]);
            double loss = population.get_loss(trip);
            if (loss<best_loss) {
                best_loss=loss;
            }
            losses1.push_back(loss);
            individui_with_losses1.emplace_back(population_1[i], loss, i);
        }

        population.sortIndividualsByLoss(individui_with_losses1);
        population.sortAscending(losses1);

        out << generation + 1 << setw(25) << calcolaMedia_first_half(losses1) << setw(25) << individui_with_losses1[0].loss << endl;
        cout << "best Loss: " << individui_with_losses1[0].loss << endl;
        individui_with_losses = individui_with_losses1;

        if (generation % 30 == 0) {
            for (int nodo = 0; nodo < size; ++nodo) {
                //if (rank == nodo) {
                    //cout << nodo << " change fatto" << endl;
                    numero_change++;
                    
                    vector<Position> exchange_trip = population.create_trip(positions, individui_with_losses[0].individual);

                    for (int i = 0; i < N_city; ++i) {
                        x_trans[i] = exchange_trip[i].get_x();
                        y_trans[i] = exchange_trip[i].get_y();
                    }
                //}

                /*// debug per il nodo che invia
                cout << "Nodo " << nodo << " sta inviando le coordinate:" << endl;
                for (int i = 0; i < 3; ++i) {
                    cout << "x[" << i << "] = " << x_trans[i] << ", y[" << i << "] = " << y_trans[i] << endl;
                }*/

                MPI_Bcast(&x_trans.front(), x_trans.size(), MPI_DOUBLE, nodo, MPI_COMM_WORLD);
                MPI_Bcast(&y_trans.front(), y_trans.size(), MPI_DOUBLE, nodo, MPI_COMM_WORLD);

                //if (rank != nodo) {
                    for (int i = 0; i < N_city; ++i) {
                        trip[i].set_x(x_trans[i]);
                        trip[i].set_y(y_trans[i]);
                    }
                //}
                // debug per i nodi che ricevono
                /*if (nodo == 0) { // Solo per il primo nodo ricevente come esempio
                    cout << "Nodo " << rank << " ha ricevuto le coordinate dal nodo " << nodo << ":" << endl;
                    for (int i = 0; i < 3; ++i) {
                        cout << "x[" << i << "] = " << x_trans[i] << ", y[" << i << "] = " << y_trans[i] << endl;
                    }
                }*/
                for (int i = 0; i < n_individui; i++) {
                    //trip = population.create_trip(positions, population_1[i]);
                    double loss = population.get_loss(trip);
                    losses1.push_back(loss);
                    individui_with_losses1.emplace_back(population_1[i], loss, i);
                }

            population.sortIndividualsByLoss(individui_with_losses1);
            population.sortAscending(losses1);

            //out << generation + 1 << setw(25) << calcolaMedia_first_half(losses1) << setw(25) << individui_with_losses1[0].loss << endl;

            individui_with_losses = individui_with_losses1;
            }
        }
    }

    ofstream outf;
    outf.open("final_trip_" + to_string(rank + 1) + ".dat");
    vector<Position> final_trip = population.create_trip(positions, individui_with_losses[0].individual);
    for (int i = 0; i < final_trip.size(); i++) {
        outf << final_trip[i].get_x() << "     " << final_trip[i].get_y() << endl;
    }
    outf.close();

    population.print_mutation_rate();
    cout << "numero change = " << numero_change << endl;

    cout << "absolute best loss = " << best_loss << endl;

    MPI_Finalize();

    return 0;
}
