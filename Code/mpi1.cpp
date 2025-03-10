#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <ostream>
#include <sstream>
#include <map>
#include <chrono>
#include <hash.hpp>
#include <geometry_basics.hpp>
#include <frechet_distance.hpp>
#define LSH_FAMILY_SIZE 8
#define LSH_SEED 234
#define LSH_RESOLUTION 80 //0.08 for taxi // generated 80
#define TAG_TASK 1
#define TAG_STOP 0
#define ODD_NODE 1
#define EVEN_NODE 2

// Parametri principali
int FILE_BUFFER_SIZE = 1000;  // Numero di linee lette e inviate per batch

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
    element_t() {

    }

    long LSH;
    int dataSet;
    curve trajectory;
    std::array<long, LSH_FAMILY_SIZE> relativeLSHs;
    long id;

    element_t(long hash,
              int dataset,
              curve trajectory,
              decltype(relativeLSHs) relativeLSHs,
              long id) : LSH(hash), dataSet(dataset), trajectory(trajectory), relativeLSHs(relativeLSHs), id(id) {}

};
//Need to create a custom type for curve to be used in MPI
MPI_Datatype createCurveType() {
    const int nitems = 3; // Numero di membri
    int blocklengths[nitems] = {TRAJ_MAX_SIZE, TRAJ_MAX_SIZE, 1}; // Dimensioni dei blocchi
    MPI_Datatype types[nitems] = {MPI_DOUBLE, MPI_DOUBLE, MPI_UNSIGNED_LONG}; // Tipi MPI
    MPI_Aint offsets[nitems];

    offsets[0] = offsetof(curve, points);         // Offset per `points`
    offsets[1] = offsetof(curve, prefix_length);  // Offset per `prefix_length`
    offsets[2] = offsetof(curve, actualSize);     // Offset per `actualSize`

    MPI_Datatype mpi_curve_type;
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_curve_type);
    MPI_Type_commit(&mpi_curve_type);

    return mpi_curve_type;
}

//Need to create a custom type for element_t to be used in MPI
MPI_Datatype createElementType() {
    MPI_Datatype mpi_curve_type = createCurveType();

    const int nitems = 5;
    int blocklengths[nitems] = {1, 1, 1, LSH_FAMILY_SIZE, 1}; // Dimensioni dei blocchi
    MPI_Datatype types[nitems] = {MPI_LONG, MPI_INT, mpi_curve_type, MPI_LONG, MPI_LONG}; // Tipi MPI
    MPI_Aint offsets[nitems];

    offsets[0] = offsetof(element_t, LSH);
    offsets[1] = offsetof(element_t, dataSet);
    offsets[2] = offsetof(element_t, trajectory);
    offsets[3] = offsetof(element_t, relativeLSHs);
    offsets[4] = offsetof(element_t, id);

    MPI_Datatype mpi_element_type;
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_element_type);
    MPI_Type_commit(&mpi_element_type);

    MPI_Type_free(&mpi_curve_type);
    return mpi_element_type;
}

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

int checkHelper(const long lsh, const element_t& a, const element_t& b){

    int foundSimilar_local = 0;

    for(size_t ii = 0; ii < a.relativeLSHs.size(); ii++)
        if (a.relativeLSHs[ii] == b.relativeLSHs[ii]){
            if (lsh == a.relativeLSHs[ii]){
                if (similarity_test(a.trajectory, b.trajectory)){
                    ++foundSimilar_local;
                    //*resultsStream << a.id << "\t" << b.id << std::endl; //TODO
                }
            }
            return foundSimilar_local;
        }

    return foundSimilar_local;
}

int main(int argc, char** argv) {

    FrechetLSH lsh_family[LSH_FAMILY_SIZE]; // 8 LSH functions
    // Initialize the LSH functions
    for (size_t i = 0; i < LSH_FAMILY_SIZE; i++)
        lsh_family[i].init(LSH_RESOLUTION, LSH_SEED*i, i); //grid_delta, seed, id
        
    MPI_Init(&argc, &argv);

    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Datatype mpi_element_type = createElementType();

    // Nodo Reader (rank 0)
    if (rank == 0) {

        std::ifstream file(argv[1]);
        if (!file.is_open()) {
            std::cerr << "Error opening dataset file!" << std::endl;
            return -1;
        }
        std::vector<std::string> fileContent;
        std::string fileline;
        while (std::getline(file, fileline)){
            //add line to vector to save content of the whole file
            fileContent.push_back(fileline);
        }
        
        //start chrono
        auto start = std::chrono::high_resolution_clock::now();

        int RoundRobin = 0;
        int totalStrings = 0;
        std::string buffer;

        MPI_Request request;
        int c = 0;
        for(auto line : fileContent){
            auto fileRead2 = std::chrono::high_resolution_clock::now();
            buffer += line + "\n";
            totalStrings++;
            if (totalStrings > FILE_BUFFER_SIZE) {
                // Invia il buffer
                int sendTo = (RoundRobin % (numProcs - 3))+3;
                MPI_Isend(buffer.data(), buffer.size(), MPI_CHAR, sendTo, TAG_TASK, MPI_COMM_WORLD, &request);
                //std::cout<<"Inviando "<<totalStrings<<" elementi a "<<(RoundRobin % (numProcs - 3))+3<<std::endl;
                RoundRobin++;
                totalStrings = 0;
                MPI_Wait(&request, MPI_STATUS_IGNORE);
                buffer.clear();
            }
        }
        MPI_Request request2;

        if (!buffer.empty()) {
            MPI_Isend(buffer.data(), buffer.size(), MPI_CHAR, (RoundRobin % (numProcs - 3))+3, TAG_TASK, MPI_COMM_WORLD, &request2);
        }
        file.close();
        //std::cout<<"File chiuso"<<std::endl;
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        // Invio segnale di fine fase 1
        for (int i = 3; i < numProcs; ++i) {
            MPI_Request request;
            MPI_Isend(nullptr, 0, MPI_CHAR, i, TAG_STOP, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
            //std::cout << "STOP to " << i << std::endl;
        }
        std::cout << "STOP EVERYONE!" << std::endl;

        //Colleziono localSimilarities dai nodi pari e dispari
        int total_similarities = 0;
        for (int i = 1; i < 3; ++i) {
            int local_similarities = 0;
            MPI_Recv(&local_similarities, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_similarities += local_similarities;
        }
        //start chrono
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Total similarities: " << total_similarities << std::endl;
        std::cout << "Elapsed time: " << std::chrono::duration<double>(end - start).count() << "s" << std::endl;
        //std::cout << "Elapsed time to read file: " << elapsed.count() << "s" << std::endl;
    }
    else if(rank == EVEN_NODE || rank == ODD_NODE){

        MPI_Status status;
        bool go = true;
        int localSimilarities = 0;
        int many = 0;
        int terminatedNodes = 0;
        std::unordered_map<long, std::vector<element_t>> localMap;

        while(go) {
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if(status.MPI_TAG == TAG_TASK) {
                //Riceve un elemento
                element_t received_element; // Buffer per ricevere l'elemento
                many++;
                MPI_Status status;
                MPI_Recv(&received_element, 1, mpi_element_type, MPI_ANY_SOURCE, TAG_TASK, MPI_COMM_WORLD, &status);
                //std::cout << "Node " << rank << "Received " << received_element.LSH << " elements" << std::endl;

                for (const element_t& insertedElement : localMap[received_element.LSH]) {
                    if (insertedElement.dataSet != received_element.dataSet) {
                        //localSimilarities += checkHelper(received_element.LSH, insertedElement, received_element);
                        int sim = checkHelper(received_element.LSH, insertedElement, received_element);
                        //if (sim > 0) {
                        //  std::cout << "Nodo " << rank << ": Similarità trovata con LSH = " << received_element.LSH << std::endl;
                        //}
                        localSimilarities += sim;

                    }
                }
                localMap[received_element.LSH].emplace_back(received_element);
            }
            //  free(tempBuffer);

            if(status.MPI_TAG == TAG_STOP) {
                MPI_Recv(nullptr, 0, MPI_CHAR, MPI_ANY_SOURCE, TAG_STOP, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //std::cout << "Node " << rank << "Received stop" << std::endl;
                terminatedNodes++;
                //Send similarities to Node 0
                if(terminatedNodes == (numProcs - 3)) {
                    MPI_Request request;
                    std::cout << "Local similarities: " << localSimilarities << " and object " << many << std::endl;
                    MPI_Isend(&localSimilarities, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
                    MPI_Wait(&request, MPI_STATUS_IGNORE);
                    /* std::cout << "Nodo " << rank << ": localMap contiene " << localMap.size() << " chiavi." << std::endl;
                     for (const auto& [key, vec] : localMap) {
                         std::cout << "LSH = " << key << ", numero di elementi = " << vec.size() << ", dataSet: ";
                         for (const auto& elem : vec) {
                             std::cout << elem.dataSet << " ";
                         }
                         std::cout << std::endl;
                     }*/
                    go = false;
                }
            }
        }
        //std::cout<<"Nodo collector terminated"<<std::endl;
    }

    else {
        std::array<std::vector<element_t>, 2> elementsToSend;
        MPI_Status status;
        bool go = true;

        while(go) {

            //std::cout<<"Node "<<rank<<" receiving.."<<std::endl;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int tag = status.MPI_TAG;
            int source = status.MPI_SOURCE;
            int count = 0;

            if (tag == TAG_TASK) {
                //std::cout << "Node " << rank << std::endl;
                // Ricezione del buffer
                MPI_Get_count(&status, MPI_CHAR, &count);
                std::vector<char> buffer(count);
                MPI_Status s;
                MPI_Recv(buffer.data(), count, MPI_CHAR, source, TAG_TASK, MPI_COMM_WORLD, &s);

                std::string buffer_str(buffer.begin(), buffer.end());
                std::istringstream stream(buffer_str);

                std::string line;
                while (std::getline(stream, line)) {
                    item parsed_item = parseLine(line);
                    //Hash computation
                    auto relative_lshs = std::array<long, LSH_FAMILY_SIZE>{};
                    for (size_t i = 0; i < LSH_FAMILY_SIZE; i++)
                        (relative_lshs)[i] = lsh_family[i].hash(parsed_item.content);

                    for (long h: relative_lshs) {
                        //std::cout<<"Hash: "<<h<<std::endl;
                        auto e = element_t(h, parsed_item.dataset, parsed_item.content, relative_lshs, parsed_item.id);
                        int pos = h % 2; //0 EVEN NODE-> 2 //1 ODD NODE->1
                        int destinationNode = pos == 0 ? EVEN_NODE : ODD_NODE;
                        //Invia elemento
                        MPI_Request request;
                        MPI_Isend(&e, 1, mpi_element_type, destinationNode, TAG_TASK, MPI_COMM_WORLD, &request);
                        MPI_Wait(&request, MPI_STATUS_IGNORE);
                    }
                }

            }

            else if (tag == TAG_STOP) {
                // Ricezione del messaggio di stop
                MPI_Recv(nullptr, 0, MPI_CHAR, source, TAG_STOP, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //Inoltro lo stop ai collectors
                MPI_Request requestEven, requestOdd;
                MPI_Isend(nullptr, 0, MPI_CHAR, EVEN_NODE, TAG_STOP, MPI_COMM_WORLD, &requestEven);
                MPI_Wait(&requestEven, MPI_STATUS_IGNORE);
                MPI_Isend(nullptr, 0, MPI_CHAR, ODD_NODE, TAG_STOP, MPI_COMM_WORLD, &requestOdd);
                MPI_Wait(&requestOdd, MPI_STATUS_IGNORE);
                go = false;
            }
            else{}
        }
        //std::cout<<"Worker terminated"<<std::endl;
    }
    //std::cout<<"Node rank "<<rank<<" Terminate"<<std::endl;
    MPI_Type_free(&mpi_element_type);
    MPI_Finalize();
    return 0;
}
