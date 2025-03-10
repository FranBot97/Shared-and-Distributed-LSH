
#include <iostream>
#include <ostream>
#include <sstream>
#include <map>
#include <vector>
#include <random>
#include <chrono>
#include "hash.hpp"
#include "geometry_basics.hpp"
#include "frechet_distance.hpp"
#include "dependencies/hash.hpp"

#define LSH_FAMILY_SIZE 8
#define LSH_SEED 234
#define LSH_RESOLUTION 80 //0.08 for taxi // generated 80

const static double SIM_THRESHOLD = 10;  //0.01 for taxi  //generated 10
const static double SIM_THRESHOLD_SQR = sqr(SIM_THRESHOLD);

size_t foundSimilar = 0;

std::ostream* resultsStream = &std::cout;

struct item {
    size_t id;
    int dataset;
    curve content;
};

struct element_t {
    long LSH;
    int dataSet;
    curve trajectory;
    std::array<long, LSH_FAMILY_SIZE> relativeLSHs;
    long id;

    element_t(long hash, int dataset, curve trajectory, decltype(relativeLSHs) relativeLSHs, long id) : LSH(hash), dataSet(dataset), trajectory(trajectory), relativeLSHs(relativeLSHs), id(id) {}
};

item parseLine(std::string& line){
    std::istringstream ss(line);
    item output;
    ss >> output.id;
    ss >> output.dataset;
    std::string tmp;
    ss >> tmp;
    if (!(tmp.find_first_of("[")==std::string::npos)) {
        tmp.replace(0, 1, "");
        tmp.replace(tmp.length() - 1, tmp.length(), "");
        bool ext = false;
        while (!ext) {
            std::string extrait = tmp.substr(tmp.find("["), tmp.find("]") + 1);
            std::string extrait1 = extrait.substr(1, extrait.find(",") - 1);
            std::string extrait2 = extrait.substr(extrait.find(",") + 1);
            extrait2 = extrait2.substr(0, extrait2.length() - 1);
            double e1 = stod(extrait1);
            double e2 = stod(extrait2);
            output.content.push_back(std::move(point(e1, e2)));
            ext = (tmp.length() == extrait.length());
            if (!ext)
                tmp = tmp.substr(tmp.find_first_of("]") + 2, tmp.length());
        }
    }
    return output;
}

bool similarity_test(const curve& c1, const curve& c2){
    // check euclidean distance
    if (euclideanSqr(c1[0], c2[0]) > SIM_THRESHOLD_SQR || euclideanSqr(c1.back(), c2.back()) > SIM_THRESHOLD_SQR)
        return false;
    // check equal time
    if (equalTime(c1, c2, SIM_THRESHOLD_SQR) || get_frechet_distance_upper_bound(c1, c2) <= SIM_THRESHOLD)
        return true;
    if (negfilter(c1, c2, SIM_THRESHOLD))
        return false;
    // full check
    if (is_frechet_distance_at_most(c1, c2, SIM_THRESHOLD))
        return true;
    return false;
}

void checkHelper(const long lsh, const element_t& a, const element_t& b){
    for(size_t ii = 0; ii < a.relativeLSHs.size(); ii++)
        if (a.relativeLSHs[ii] == b.relativeLSHs[ii]){
            if (lsh == a.relativeLSHs[ii]){
                if (similarity_test(a.trajectory, b.trajectory)){
                    ++foundSimilar;
                    //*resultsStream << a.id << "\t" << b.id << std::endl;
                }
            }
            return;
        }
}

int main(int argc, char** argv){

    if (argc < 2){
        std::cout << "Usage: " << argv[0] << " inputDataset <outputFile>";
        return EXIT_FAILURE;
    }

    std::ofstream filestream;
    if (argc > 2){
        filestream = std::ofstream(argv[2]);
        if (filestream.is_open())
            resultsStream = &filestream;
    }
    std::string line;

    std::ifstream file(argv[1]);

    // Check if the file is opened successfully
    if (!file.is_open()) {
        std::cerr << "Error opening dataset file!" << std::endl;
        return -1;
    }

    // build the LSH function family
    // Vengono inizializzate 8 diverse funzioni LSH, ognuna differisce per il seed (e id)
    // nella funzione this->hash viene calcolato il valore della funzione LSH per la curva
    // usando i parametri inizializzati per ottenere il un valore di hash
    FrechetLSH lsh_family[LSH_FAMILY_SIZE]; // 8 LSH functions
    // Initialize the LSH functions
    for (size_t i = 0; i < LSH_FAMILY_SIZE; i++)
        lsh_family[i].init(LSH_RESOLUTION, LSH_SEED*i, i); //grid_delta, seed, id

    std::unordered_map<long, std::vector<element_t>> elements;
    std::chrono::duration<double> elapsed{};
    std::vector<string> lines;

    while(std::getline(file,line)){
        lines.emplace_back(line);
    }
    // Close the file
    file.close();

    //Get time
    auto start = std::chrono::high_resolution_clock::now();

    for(auto l : lines){
        // Process each line here
        item it = parseLine(l);
        // Compute the LSH values for the curve
        std::array<long, LSH_FAMILY_SIZE> relative_lshs;
        for(size_t i = 0; i < LSH_FAMILY_SIZE; i++)
            relative_lshs[i] = lsh_family[i].hash(it.content);

        for(long h : relative_lshs)
            elements[h].emplace_back(h, it.dataset, it.content, relative_lshs, it.id);

    }

    for (auto& [lsh, elements_v] : elements)
        for(size_t i = 0; i < elements_v.size(); i++)
            for(size_t j = i+1; j < elements_v.size(); j++)
                if (elements_v[i].dataSet != elements_v[j].dataSet)
                    checkHelper(lsh, elements_v[i],elements_v[j]);

    //Calculate time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> totalElapsed = end - start;

    std::cout << "Trajectories similar found: " << foundSimilar << std::endl;
    std::cout << "Total time after read: " << totalElapsed.count() << "s" << std::endl;

    return 0;
}
