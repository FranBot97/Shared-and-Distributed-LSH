#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <ostream>
#include <sstream>
#include <cstring>
#include <map>
#include <chrono>
#include <dependencies/hash.hpp>
#include <dependencies/geometry_basics.hpp>
#include <dependencies/frechet_distance.hpp>
#include <cassert>
#include <omp.h>


#define LSH_FAMILY_SIZE 8
#define LSH_SEED 234
#define LSH_RESOLUTION 80 //0.08 for taxi // generated 80
#define BATCH_SIZE 1000
const static double SIM_THRESHOLD = 10;  //0.01 for taxi  //generated 10
const static double SIM_THRESHOLD_SQR = sqr(SIM_THRESHOLD);

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

    element_t() = default;

    element_t(long hash,
              int dataset,
              curve trajectory,
              decltype(relativeLSHs) relativeLSHs,
              long id) : LSH(hash), dataSet(dataset), trajectory(trajectory), relativeLSHs(relativeLSHs), id(id) {}

};

void printElement(const element_t& e, int rank) {
    std::cout << "Processo " << rank
              << " | ID: " << e.id
              << " | LSH: " << e.LSH
              << " | Dataset: " << e.dataSet
              << " | Hashes: [";
    for (size_t i = 0; i < e.relativeLSHs.size(); ++i) {
        std::cout << e.relativeLSHs[i];
        if (i < e.relativeLSHs.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << "Processo " << rank << " | Trajectory: " << e.trajectory << std::endl;
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
    }return output;
}

bool similarity_test(const curve& c1, const curve& c2){
    if (euclideanSqr(c1[0], c2[0]) > SIM_THRESHOLD_SQR || euclideanSqr(c1.back(), c2.back()) > SIM_THRESHOLD_SQR)
        return false;
    if (equalTime(c1, c2, SIM_THRESHOLD_SQR) || get_frechet_distance_upper_bound(c1, c2) <= SIM_THRESHOLD)
        return true;
    if (negfilter(c1, c2, SIM_THRESHOLD))
        return false;
    if (is_frechet_distance_at_most(c1, c2, SIM_THRESHOLD))
        return true;
    return false;
}

int checkHelper(const long lsh, const element_t& a, const element_t& b, int rank){
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

std::vector<char> serialize_element(const element_t& elem) {
    std::vector<char> buffer(sizeof(element_t));
    std::memcpy(buffer.data(), &elem, sizeof(element_t));
    return buffer;
}

element_t deserialize_element(const std::vector<char>& buffer) {
    element_t elem;
    std::memcpy(&elem, buffer.data(), sizeof(element_t));
    return elem;
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
/*    MPI_Datatype mpi_element_type = createElementType();*/

    std::vector<int> send_counts(numProcs, 0);  // Numero di caratteri da inviare a ciascun processo
    std::vector<int> displacements(numProcs, 0);
    std::vector<int> line_sizes;  // Per tracciare la lunghezza di ogni stringa
    std::vector<char> send_buffer;  // Buffer unico per tutte le stringhe concatenate
    std::vector<char> recv_buffer;  // Buffer per la ricezione

    // **Nodo 0 legge il file e prepara il buffer**
    int max_iter = 0;
    std::unordered_map<long, std::vector<element_t>> localMap;
    int localSimilarities = 0;
    int totalReceived = 0;
    std::chrono::high_resolution_clock::time_point start;
    std::vector<std::string> lines;
    std::string line;

    if (rank == 0) {
        std::ifstream file(argv[1]);
        if (!file.is_open()) {
            std::cerr << "Errore nell'apertura del file!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        std::cout << "Reading file.." << std::endl;
        while (std::getline(file, line))
            lines.push_back(line + "\n");  // Aggiungiamo `\n` per mantenere il formato
        file.close();
        std::cout << "Completed reading file" << std::endl;
    }
        start = std::chrono::high_resolution_clock::now();

        while (true) {
            auto computation = std::chrono::high_resolution_clock::now();
            if (rank == 0) {
                std::vector<std::string> partial_lines;
                max_iter = 0;
                size_t n = BATCH_SIZE * numProcs; // Numero di righe da inviare per processo
                size_t num_to_move = std::min(n, lines.size());
                partial_lines.reserve(num_to_move);
                for (size_t i = 0; i < num_to_move; ++i)
                    partial_lines.push_back(std::move(lines[i])); // Spostiamo le stringhe
                // Rimuoviamo le righe già elaborate
                lines.erase(lines.begin(), lines.begin() + num_to_move);
                int num_lines = partial_lines.size();
                //std::cout << "Process " << rank << " sending " << num_lines << " lines" << std::endl;
                int base_lines = num_lines / numProcs;
                int extra_lines = num_lines % numProcs;

                std::vector<std::string> send_strings(numProcs);
                int current_index = 0;
                for (int i = 0; i < numProcs; i++) {
                    int lines_for_this_proc = base_lines + (i < extra_lines ? 1 : 0);
                    if (lines_for_this_proc > max_iter)
                        max_iter = lines_for_this_proc;
                    for (int j = 0; j < lines_for_this_proc; j++) {
                        send_strings[i] += partial_lines[current_index++];
                    }
                    send_counts[i] = send_strings[i].size();  // Numero di caratteri per ogni processo
                }

                int total_size = 0;
                for (int count: send_counts)
                    total_size += count;
                send_buffer.resize(total_size);

                // **Copiamo le stringhe nel buffer unico senza `memcpy`**/
                int offset = 0;
                for (int i = 0; i < numProcs; i++) {
                    std::memcpy(&send_buffer[offset], send_strings[i].c_str(), send_counts[i]);
                    displacements[i] = offset;
                    offset += send_counts[i];
                }
                //std::cout<<"Process "<<rank<<" max_iter"<<max_iter<<std::endl;
            }
            //Inviamo il numero max di iterazioni
            MPI_Bcast(&max_iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (max_iter == 0)
                break;
            // **Distribuiamo la dimensione dei dati che ogni processo riceverà**
            int recv_count;
            MPI_Scatter(send_counts.data(), 1, MPI_INT, &recv_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
            recv_buffer.resize(recv_count);
            MPI_Scatterv(send_buffer.data(), send_counts.data(), displacements.data(), MPI_CHAR,
                         recv_buffer.data(), recv_count, MPI_CHAR, 0, MPI_COMM_WORLD);
            std::stringstream ss(std::string(recv_buffer.begin(), recv_buffer.end()));
            std::string received_line;
            std::vector<std::string> received_lines;
            while (std::getline(ss, received_line)) {
                received_lines.push_back(received_line);
            }
            //std::cout << "Process " << rank << " received " << received_lines.size() << std::endl;

            std::vector<element_t> local_elements;
            std::vector<std::vector<element_t>> send_elements(numProcs);
            int count = 0;

#pragma omp parallel for
            for (int i = 0; i < max_iter; i++) {
                if (i < received_lines.size()) {
                    std::string line = received_lines[i];
                    item it = parseLine(line);
                    std::array<long, LSH_FAMILY_SIZE> relative_lshs{};
                    for (size_t j = 0; j < LSH_FAMILY_SIZE; j++)
                        relative_lshs[j] = lsh_family[j].hash(it.content);
                    for (long h: relative_lshs) {
                        auto e = element_t(h, it.dataset, it.content, relative_lshs, it.id);
                        long pos = h % numProcs;
                        if (pos == rank) {
#pragma omp critical
                            local_elements.emplace_back(e);
                        } else {
#pragma omp critical
                            send_elements[pos].emplace_back(e);
                        }
                    }
                }
            }
            

            // **Alltoall dopo max_iter lines**//
            std::vector<int> send_counts_elem(numProcs, 0);
            std::vector<int> send_displacements_elem(numProcs, 0);
            std::vector<char> send_buffer_elem;
            // Serializzazione e creazione buffer di invio
            for (int j = 0; j < numProcs; j++) {
                for (const auto &elem: send_elements[j]) {
                    auto serialized = serialize_element(elem);
                    send_buffer_elem.insert(send_buffer_elem.end(), serialized.begin(), serialized.end());
                    send_counts_elem[j] += serialized.size();
                }
            }
            //std::cout<< "Process " << rank << " is sending TOTAL " << send_buffer_elem.size() << " elements\n";
            send_displacements_elem[0] = 0; // First displacement is always 0
            for (int j = 1; j < numProcs; j++)
                send_displacements_elem[j] = send_displacements_elem[j - 1] + send_counts_elem[j - 1];
            std::vector<int> recv_counts_elem(numProcs, 0);
            //std::cout << "All to all " << rank << std::endl;
            auto computation_end = std::chrono::high_resolution_clock::now();
        //    std::cout<<"Process "<<rank<<" computation time "<<std::chrono::duration_cast<std::chrono::milliseconds>(computation_end - computation).count()<<std::endl;
            MPI_Alltoall(send_counts_elem.data(), 1, MPI_INT, recv_counts_elem.data(), 1, MPI_INT, MPI_COMM_WORLD);
            int total_recv = 0;
            for (int j = 0; j < numProcs; j++)
                total_recv += recv_counts_elem[j];
            //std::cout << "Process " << rank << " is receiving " << total_recv << " elements\n";
            std::vector<char> recv_buffer_elem(total_recv);
            totalReceived += total_recv;
            std::vector<int> recv_displacements_elem(numProcs, 0);
            recv_displacements_elem[0] = 0; // First displacement is always 0
            for (int j = 1; j < numProcs; j++)
                recv_displacements_elem[j] = recv_displacements_elem[j - 1] + recv_counts_elem[j - 1];
            MPI_Request request;
            MPI_Ialltoallv(send_buffer_elem.data(), send_counts_elem.data(), send_displacements_elem.data(), MPI_CHAR,
                           recv_buffer_elem.data(), recv_counts_elem.data(), recv_displacements_elem.data(), MPI_CHAR,
                           MPI_COMM_WORLD, &request);
            for (auto e: local_elements) {
                assert(e.LSH % numProcs == rank);
                for (auto other_element: localMap[e.LSH]) {
                    if (e.dataSet != other_element.dataSet) {
                        int sim = checkHelper(e.LSH, other_element, e, rank);
                        localSimilarities += sim;
                    }
                }
                localMap[e.LSH].emplace_back(e);
            }
            local_elements.clear();
            MPI_Wait(&request, MPI_STATUS_IGNORE);
            std::vector<element_t> received_elements;
            size_t offset = 0;
            while (offset < recv_buffer_elem.size()) {
                std::vector<char> chunk(recv_buffer_elem.begin() + offset,
                                        recv_buffer_elem.begin() + offset + sizeof(element_t));
                received_elements.push_back(deserialize_element(chunk));
                offset += sizeof(element_t);
            }
            for (const auto &e: received_elements) {
                assert(e.LSH % numProcs == rank);
                for (auto other_element: localMap[e.LSH]) {
                    if (e.dataSet != other_element.dataSet) {
                        int sim = checkHelper(e.LSH, other_element, e, rank);
                        localSimilarities += sim;
                    }
                }
                localMap[e.LSH].emplace_back(e);
            }
            for (int j = 0; j < numProcs; j++)
                send_elements[j].clear();
        received_lines.clear();
    }
    std::vector<int> allSimilarities; // Root process will store all similarities here
    if (rank == 0)
        allSimilarities.resize(numProcs); // Root process allocates space for all similarities
    MPI_Gather(&localSimilarities, 1, MPI_INT, allSimilarities.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    int totalSimilarities = 0;
    if (rank == 0) {
        // Print or process the gathered similarities
        for (int i = 0; i < numProcs; i++) {
           totalSimilarities += allSimilarities[i];
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Total Similarities: "<< totalSimilarities << " similarities." << std::endl;
        std::chrono::duration<double> totalElapsed = end - start;
        std::cout << "Total time: " << totalElapsed.count() << "s" << std::endl;
    }
    MPI_Finalize();
    return 0;
}
