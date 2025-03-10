
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
#include <ff/ff.hpp>
#include <ff/parallel_for.hpp>
#include <ff/map.hpp>


#define LSH_FAMILY_SIZE 8
#define LSH_SEED 234
#define LSH_RESOLUTION 80 //0.08 for taxi // generated 80
#define BUFFER_SIZE 1000

const static double SIM_THRESHOLD = 10; //0.01;  //0.01 for taxi  //generated 10
const static double SIM_THRESHOLD_SQR = sqr(SIM_THRESHOLD);

//GLOBAL VARIABLES
int globalSimilarities = 0;
//int nWorkers = 0;
int nWorkers_Farm1 = 0;
int nWorkers_Farm2 = 0;
FrechetLSH lsh_family[LSH_FAMILY_SIZE];
std::ostream* resultsStream = &std::cout;
uint* WorkersID = nullptr;
std::chrono::time_point<std::chrono::high_resolution_clock> globalTimeStart;
std::vector<string> lines;

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

    element_t(long hash,
              int dataset,
              curve trajectory,
              decltype(relativeLSHs) relativeLSHs,
              long id) :
            LSH(hash),
            dataSet(dataset),
            trajectory(trajectory),
            relativeLSHs(relativeLSHs),
            id(id) {}
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

int checkHelper(const long lsh, const element_t& a, const element_t& b){

    int foundSimilar_local = 0;

    for(size_t ii = 0; ii < a.relativeLSHs.size(); ii++)
        if (a.relativeLSHs[ii] == b.relativeLSHs[ii]){
            if (lsh == a.relativeLSHs[ii]){
                if (similarity_test(a.trajectory, b.trajectory)){
                    ++foundSimilar_local;
                    //*resultsStream << a.id << "\t" << b.id << std::endl;
                }
            }
            return foundSimilar_local;
        }

    return foundSimilar_local;
}

/******** FARM 1 ********/
struct Dealer : ff::ff_monode_t<std::vector<std::string>> {

    std::vector<std::string>* svc(std::vector<std::string>*) override {

        auto* buffer = new std::vector<std::string>;
        buffer->reserve(BUFFER_SIZE);

        for(const auto& l : lines){
            buffer->emplace_back(l);
            if (buffer->size() >= BUFFER_SIZE) {
                ff_send_out(buffer); // Send the buffer to workers
                buffer = new std::vector<std::string>;
                buffer->reserve(BUFFER_SIZE);
            }
        }
        if (!buffer->empty()) {
            ff_send_out(buffer);
        } else {
            delete buffer;
        }
        return EOS; // Signal end of stream
    }

};

struct Parser: ff::ff_node_t<std::vector<std::string>>{
    int id;
    explicit Parser(int id): id(id) {

    }
    std::chrono::duration<double> Worker_time = std::chrono::duration<double>(0);
    std::vector<item> *items;

    vector<std::string>* svc(vector<std::string>* t) override {
        items = new std::vector<item>();
        //start timer
        auto start = std::chrono::high_resolution_clock::now();
        // Cast the input to a buffer of lines
        auto* lines = static_cast<std::vector<std::string>*>(t);

        // Process each line in the buffer
        for (auto& line : *lines) {
            // Parse the line into an item
            item it = item(parseLine(line));
            items->push_back(it);
        }
        ff_send_out(items);

        auto end = std::chrono::high_resolution_clock::now();
        Worker_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        //std::cout<<"Time to compute all lines"<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

        // Clean up the buffer
        delete lines;
        return GO_ON;
    }

    void svc_end() override {
        // std::cout<<"Time for worker:"<<Worker_time.count()<<std::endl;
        //take time
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - globalTimeStart;
        //std::cout<<"Farm1 duration: "<<elapsed.count()<<"s"<<std::endl;
    }

};

struct Router: ff::ff_monode_t<std::vector<item>>{

    std::chrono::duration<double> Worker_time = std::chrono::duration<double>(0);

    std::vector<item>* svc(std::vector<item>* t) override {
        auto *items = static_cast<std::vector<item> *>(t);
        // {   std::lock_guard<std::mutex> lock(mtx_elements);
        for (auto it: *items) {
            auto relative_lshs = std::array<long, LSH_FAMILY_SIZE>{};
            for (size_t i = 0; i < LSH_FAMILY_SIZE; i++)
                (relative_lshs)[i] = lsh_family[i].hash(it.content);

            // Crea il nuovo elemento a partire dall'oggetto item
            {   //std::lock_guard<std::mutex> lock(mtx_elements);
                for (long h: relative_lshs) {
                    auto *e = new element_t(h, it.dataset, it.content, relative_lshs, it.id);
                    //in base al valore della hash router manda l'elemento al worker corrispondente
                    ff_send_out_to(e, WorkersID[h % nWorkers_Farm2]);
                }
            }
        }

        delete items;

        return GO_ON;
    }
};


struct Counter: ff::ff_node_t<element_t> {

    long calculated = 0;
    int localSimilarities = 0;
    uint position;
    //local unordered map
    std::unordered_map<long, std::vector<element_t*>> localMap;
    //costruttore inizializza position
    explicit Counter(ulong position) : position(position) {}

    int svc_init() override {
        WorkersID[position] = get_my_id();
        return 0;
    }

    element_t* svc (element_t* e) override {
        calculated++;
        //vettore corrispondente
        //confronto l'elemento appena arrivato con tutti quelli che ci sono già
        for (const element_t* insertedElement: localMap[e->LSH]) {

            if (insertedElement->dataSet != e->dataSet) {
                localSimilarities += checkHelper(e->LSH, *insertedElement, *e);
            }
        }
        localMap[e->LSH].emplace_back(e);
        return GO_ON;
    }

/*    void svc_end() override {
        for (auto& [key, local_elements] : localMap) {
            for (element_t* elem : local_elements) {
                delete elem;
            }
        }
        localMap.clear();
    }*/

    void eosnotify(ssize_t) override {
        //std::cout<<"Im worker "<<get_my_id()<<" and I elaborated "<<calculated<<" similarities"<<std::endl;
        ff_send_out(new int(localSimilarities));
    }

};

struct FinalCollector : ff::ff_minode_t<int> {

    int WorkersFinished = 0;

    int* svc(int* t) override{
        globalSimilarities += *t;
        delete t;

        if(++WorkersFinished == nWorkers_Farm2)
            return EOS;

        return GO_ON;
    }

/*    void svc_end() override {
        *resultsStream << globalSimilarities << std::endl;
    }*/
};


int main(int argc, char** argv){

    if (argc < 4){
        std::cout << "Usage: " << argv[0] << " inputDataset nWorkersFarm1 nWorkersFarm2 <outputFile>";
        return EXIT_FAILURE;
    }

    std::ofstream filestream;
    if (argc > 4){
        filestream = std::ofstream(argv[4]);
        if (filestream.is_open())
            resultsStream = &filestream;
    }

    nWorkers_Farm1 = std::stoi(argv[2]);
    nWorkers_Farm2 = std::stoi(argv[3]);
    WorkersID = new uint[nWorkers_Farm2];

    std::ifstream file(argv[1]);
    if (!file.is_open()) {
        std::cerr << "Error opening dataset file!" << std::endl;
        return -1;
    }

    std::string line;

    while(std::getline(file,line)){
        lines.emplace_back(line);
    }

    //START COUNTING
    std::cout << "Start timer.." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    //We can avoid parallel for if LSH_FAMILY_SIZE is small
    ff::ParallelFor pf(nWorkers_Farm1);
    pf.parallel_for(0, LSH_FAMILY_SIZE, 1, [&](const long i){
        lsh_family[i].init(LSH_RESOLUTION, LSH_SEED*i, i);
    });

    /** FARM 1 - STAGE 1 **/
    //FileReader -> Parser[0..nWorkers]
    Dealer LinesManager;
    std::vector<std::unique_ptr<ff::ff_node>> Workers_Stage1;
    for (int i = 0; i < nWorkers_Farm1; ++i) {
        Workers_Stage1.push_back(make_unique<Parser>(i));
    }
    ff::ff_Farm<> Farm_Stage1(std::move(Workers_Stage1));
    Farm_Stage1.add_emitter(LinesManager);
    Farm_Stage1.remove_collector();

    /** FARM 2 - STAGE 2 **/
    //Router -> Counter[0..nWorkers]
    Router Router;
    std::vector<std::unique_ptr<ff::ff_node>> Workers_Stage2;
    for (int i = 0; i < nWorkers_Farm2; ++i) {
        Workers_Stage2.push_back(make_unique<Counter>(i));
    }
    ff::ff_Farm<> Farm_Stage2(std::move(Workers_Stage2));
    FinalCollector FinalCollector;
    Farm_Stage2.add_emitter(Router);
    Farm_Stage2.add_collector(FinalCollector);

    ff::ff_Pipe LSHSimilarityPipe(Farm_Stage1, Farm_Stage2);

    globalTimeStart = std::chrono::high_resolution_clock::now();
    if (LSHSimilarityPipe.run_and_wait_end()<0) {
        std::cout<<"error running pipe";
        return -1;
    }

    //Calculate time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> totalElapsed = end - start;

    std::cout << "Trajectories similar found: " << globalSimilarities << std::endl;
    std::cout << "Total time after read: " << totalElapsed.count() << "s" << std::endl;


    return 0;
}
