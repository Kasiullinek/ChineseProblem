#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

struct switches {
    std::string i{};
    std::string o{};
    int p{};
};

struct inputFile {
    int intersection1{};
    int intersection2{};
    double streetLength{};
    std::string streetName{};
};

typedef std::map <int, std::set<std::pair<int, double>>> Graph;

void SaveNew(std::string& x, const std::string& fileToSave) {
    std::ofstream out{ fileToSave };
    if (out) {
        out << x << std::endl;
    }
    out.close();
}

std::string LoadNew(const std::string& fileWithSavedX) {
    std::string x;
    std::ifstream in(fileWithSavedX);
    if (in) {
        std::string line;
        while (std::getline(in, line)) {
            if (line.length() == 0) break;
            std::stringstream ss(line);
            std::string newFileName;
            std::getline(ss, newFileName, ',');
            x = newFileName;
        }
        in.close();
    }
    return x;
}

std::vector<inputFile> LoadFromFile(std::string& fileName, bool& isFileLoaded, const std::string& fileWithDecision) {

    std::vector<inputFile> result;
    std::ifstream in(fileName);

    if (in) {
        isFileLoaded = true;
        std::string line;

        while (std::getline(in, line)) {

            if (line.length() == 0) break;

            std::stringstream ss(line);
            std::string inter1, inter2, length, name;

            std::getline(ss, inter1, ' ');
            std::getline(ss, inter2, ' ');
            std::getline(ss, length, ' ');
            std::getline(ss, name, ',');

            result.push_back({ stoi(inter1), stoi(inter2), stod(length), name });
        }
        in.close();
    }
    else {
        isFileLoaded = false;
        std::cout << std::endl;
        std::cout << "Nie znaleziono pliku wejsciowego.";
    }

    std::string stringIsFileLoaded = std::to_string(isFileLoaded);
    SaveNew(stringIsFileLoaded, fileWithDecision);

    return result;
}

void PrintInputFile(const std::vector<inputFile>& vector) {
    std::cout << std::endl;
    for (const auto& el : vector) {
        std::cout << el.intersection1 << " " << el.intersection2 << " " << el.streetLength << " " << el.streetName << std::endl;
    }
    std::cout << std::endl;
}

bool IsStartPoint(const int point, const std::vector<inputFile>& vector) {
    for (const auto& el : vector) {
        if ((point == el.intersection1) || (point == el.intersection2)) {
            return true;
        }
    }
    return false;
}

Graph LoadGraph(const std::vector<inputFile>& vector) {
    Graph graph;
    for (const auto& el : vector) {
        graph[el.intersection1].insert({ el.intersection2, el.streetLength });
        graph[el.intersection2].insert({ el.intersection1, el.streetLength });
    }
    return graph;
}

void AddEdge(Graph& graph, const int& node1, const int&node2, const int& length) {
    graph[node1].insert({node2, length});
    graph[node2].insert({node1, length});
}

//robocze
void PrintGraph(const Graph& graph) {

    for (auto i = graph.begin(); i != graph.end(); i++) {
        auto nodeI = i->first;
        std::cout << "Node " << nodeI << " makes an edge with " << std::endl;
        std::set<std::pair<int, double>> nodeJ_Length = i->second;
        for (const auto& el : nodeJ_Length) {
            std::cout << "\tNode " << el.first << " with edge weight =" << el.second << std::endl;
        }
        std::cout << std::endl;
    }
}

std::vector<std::vector<int>> LoadMatrix(Graph& graph) {
    size_t graphSize = graph.size();

    std::vector< std::vector <int> > matrix;
    for (size_t i = 0; i < graphSize; i++) {
        std::vector <int> row;
        for (size_t j = 0; j < graphSize; j++) {
            row.push_back(0);
        }
        matrix.push_back(row);
    }
    for (size_t i = 0; i < graphSize; i++) {
        for (size_t j = 0; j < graphSize; j++) {
            matrix[i][j] = 0;
        }
    }
    for (auto i = graph.begin(); i != graph.end(); i++) {
        int nodeE = i->first;
        int e = (nodeE)-1;
        std::set<std::pair<int, double>> nodeE_Length = i->second;
        for (const auto& el : nodeE_Length) {
            int nodeF = el.first;
            int f = (nodeF)-1;
            matrix[e][f]++;
            matrix[f][e]++;
        }
    }
    return matrix;
}

std::vector<std::vector<int>> LoadMatrix2(const std::vector<inputFile>& vector,const Graph& graph) {
    size_t vectorSize = vector.size();
    size_t graphSize = graph.size();

    std::vector< std::vector <int> > matrix;
    for (size_t i = 0; i < graphSize; i++) {
        std::vector <int> row;
        for (size_t j = 0; j < graphSize; j++) {
            row.push_back(0);
        }
        matrix.push_back(row);
    }
    for (size_t i = 0; i < graphSize; i++) {
        for (size_t j = 0; j < graphSize; j++) {
            matrix[i][j] = 0;
        }
    }
    for (auto i = graph.begin(); i != graph.end(); i++) {
        int nodeE = i->first;
        int e = (nodeE)-1;
        std::set<std::pair<int, double>> nodeE_Length = i->second;
        for (const auto& el : nodeE_Length) {
            int nodeF = el.first;
            int streetLength = int(el.second);
            int f = (nodeF)-1;
            matrix[e][f] = streetLength;
            matrix[f][e] = streetLength;
        }
    }
    return matrix;
}

//robocze
void PrintMatrix(const std::vector<std::vector<int>>& matrix, const int& graphSize) {
    for (int i = 0; i < graphSize; i++) {
        for (int j = 0; j < graphSize; j++) {
            std::cout << "|" << matrix[i][j] << "|";
        }
        std::cout << std::endl;
    }
}

void Traverse(size_t& n, bool visited[], const size_t& numberOfNodes, const size_t& numberOfStreets, const std::vector<std::vector<int>>& matrix) {

    visited[n] = true;

    for (size_t i = 0; i < numberOfNodes; i++) {
        if (matrix[n][i]) {
            if (!visited[i])
                Traverse(i, visited, numberOfNodes, numberOfStreets, matrix);
        }
    }
}

bool IsGraphConnected(const size_t& numberOfStreets, const size_t& numberOfNodes, const std::vector<std::vector<int>>& matrix) {
    bool* visited = new bool[numberOfNodes];

    for (size_t i = 0; i < numberOfNodes; i++) {
        for (size_t j = 0; j < numberOfNodes; j++) {
            visited[j] = false;
        }

        Traverse(i, visited, numberOfNodes, numberOfStreets, matrix);

        for (size_t i = 0; i < numberOfNodes; i++) {
            if (!visited[i]) {
                return false;
            }
        }
    }
    return true;
}

int CountOddNodes(const Graph& graph) {
    size_t numberOfNodes = graph.size();
    std::vector<int> listOfNodes(numberOfNodes + 1);
    int oddNodesCounter = 0;

    for (auto i = graph.begin(); i != graph.end(); i++) {
        auto nodeI = i->first;
        std::set<std::pair<int, double>> nodeJ_Length = i->second;
        size_t numberOfNodes = nodeJ_Length.size();
        listOfNodes[nodeI] = int(numberOfNodes);
        if ((listOfNodes[nodeI] % 2) != 0) {
            oddNodesCounter++;
        }
    }
    return oddNodesCounter;
}

std::vector<int> LoadOddNodes(const Graph& graph, const int& numberOfOddNodes) {
    size_t numberOfNodes = graph.size();
    std::vector<int> listOfNodes(numberOfNodes + 1);
    std::vector<int> listOfOddNodes(numberOfOddNodes);
    int oddNodesCounter = 0;

    for (auto i = graph.begin(); i != graph.end(); i++) {
        auto nodeI = i->first;
        std::set<std::pair<int, double>> nodeJ_Length = i->second;
        size_t numberOfNodes = nodeJ_Length.size();
        listOfNodes[nodeI] = int(numberOfNodes);
        if ((listOfNodes[nodeI] % 2) != 0) {
            listOfOddNodes[oddNodesCounter] = nodeI;
            oddNodesCounter++;
        }
    }
    return listOfOddNodes;
}

//robocze
void PrintEulerCircle(const std::vector<std::pair<int, int>>& list) {
    for (const auto& el : list) {
        std::cout << el.first << "-" << el.second << " ";
    }
}

std::vector<std::pair<int, int>> LoadEulerCycle(const std::string& fileWithEulerCycle) {
    std::vector<std::pair<int, int>> result;
    std::ifstream in(fileWithEulerCycle);

    if (in) {
        std::string line;

        while (std::getline(in, line)) {

            if (line.length() == 0) break;

            std::stringstream ss(line);
            std::string inter1, inter2;

            std::getline(ss, inter1, ' ');
            std::getline(ss, inter2, ' ');

            result.push_back({ stoi(inter1), stoi(inter2) });
        }
        in.close();
    }
    return result;
}

void SaveEulerCycle(int& node1, int& node2, std::ofstream& out) {
    out << node1 << " " << node2 << std::endl;
}

//
bool IsBridge(int u, int v, const size_t& numberOfNodes, std::vector<std::vector<int>>& matrix) {
    int degree = 0;
    for (size_t i = 0; i < numberOfNodes; i++) {
        if (matrix[v][i]) {
            degree++;
        }
    }
    if (degree > 1) {
        return false;
    }
    return true;
}

//
int EdgeCount(std::vector<std::vector<int>>& matrix, const size_t& numberOfNodes) {
    int count = 0;
    for (size_t i = 0; i < numberOfNodes; i++)
        for (size_t j = i; j < numberOfNodes; j++)
            if (matrix[i][j])
                count++;
    return count;
}

//
void FleuryAlgorithm(int& startPoint, size_t& numberOfStreets, const size_t& numberOfNodes, std::vector<std::vector<int>>& matrix, std::ofstream& out) {
    static int edge = EdgeCount(matrix, numberOfNodes);
    for (size_t i = 0; i < numberOfNodes; i++) {
        if (matrix[startPoint][i]) {
            if (edge <= 1 || !IsBridge(startPoint, i, numberOfNodes, matrix)) {
                int node1 = startPoint + 1;
                int node2 = i + 1;
                SaveEulerCycle(node1, node2, out);
                matrix[startPoint][i] = matrix[i][startPoint] = 0;
                edge--;
                int x = int(i);
                FleuryAlgorithm(x, numberOfStreets, numberOfNodes, matrix, out);
            }
        }
    }
}

void SaveToOutputFile(int& intersection1, int& intersection2, std::string& streetName, std::ofstream& fout) {
    fout << intersection1 << " " << intersection2 << " " << streetName << std::endl;
}

//
void MakeOutputFile(std::vector<std::pair<int, int>>& e, std::vector<inputFile>& v, int& iterator, std::ofstream& fout) {
    size_t numberOfNodes = v.size();
    size_t numberOfPairs = e.size() - 1;
    for (size_t i = 0; i < numberOfNodes; i++) {
        if (((e[iterator].first == v[i].intersection1) &&
            (e[iterator].second == v[i].intersection2)) ||
            ((e[iterator].first == v[i].intersection2) &&
                (e[iterator].second == v[i].intersection1))) {
            std::cout << e[iterator].first << " " << e[iterator].second << " " << v[i].streetName << std::endl;
            SaveToOutputFile(e[iterator].first, e[iterator].second, v[i].streetName, fout);
        }
    }

    iterator++;
    if (iterator <= int(numberOfPairs)) {
        MakeOutputFile(e, v, iterator, fout);
    }
    else {
        return;
    }
}

std::vector<std::pair<int, int>> CombinePairs(const std::vector<int> listOfOddNodes) {
    int listSize = listOfOddNodes.size();
    int nrOfCombinations = (listSize*(listSize-1))/2;
    int* listOfOddNodes2 = new int [listSize];
    std::vector<std::pair<int, int>> listOfPairs(nrOfCombinations);
    int iterator = 0;

    for (int i = 0; i < listSize; i++) {
        listOfOddNodes2[i] = listOfOddNodes[i];
    }

    for (int i = 0; i < listSize; i++) {
        for (int j = 0; j < listSize; j++) {
            if (i != j && i < j) {
                listOfPairs[iterator].first = listOfOddNodes[i];
                std::cout << listOfPairs[iterator].first << " ";
                listOfPairs[iterator].second = listOfOddNodes2[j];
                std::cout << listOfPairs[iterator].second << std::endl;
                iterator++;
            }
        }
    }

    return listOfPairs;
    
    delete[] listOfOddNodes2;
}

void DijkstraAlgorithm(const std::vector<std::vector<int>>& matrix2, const int& startNode, const int& endNode) {
    int matrix2Size = matrix2.size();
    int* distance = new int[matrix2Size];
    int* pred = new int[matrix2Size];
    int* pred2 = new int[matrix2Size];
    bool* visited = new bool[matrix2Size];

    for (int i = 0; i < matrix2Size; i++) {
        distance[i] = INT_MAX;
        pred[i] = startNode;
        pred2[i] = startNode;
        visited[i] = false;
    }
    distance[startNode] = 0;

    for (int count = 0; count < matrix2Size - 1; count++) {

        int minDistance = INT_MAX;
        int nextNode;
        for (int i = 0; i < matrix2Size; i++) {
            if (visited[i] == false && distance[i] <= minDistance) {
                minDistance = distance[i];
                nextNode = i;
            }
        }
        visited[nextNode] = true;
        for (int i = 0; i < matrix2Size; i++) {
            if (!visited[i] && matrix2[nextNode][i] && distance[nextNode] != INT_MAX &&
                distance[nextNode] + matrix2[nextNode][i] < distance[i]) {
                distance[i] = distance[nextNode] + matrix2[nextNode][i];
                pred[i] = nextNode;
                pred2[i] = nextNode;
            }

        }
    };

    std::pair<std::vector<int>, std::vector<double>> sp;

    std::cout << "\nDistance from " << startNode + 1 << " to node " << endNode + 1 << "=" << distance[endNode];
    std::cout << "\nPath=" << endNode + 1;
    int j = endNode;
    int k = endNode;
    do {
        j = pred[j];
        std::cout << "<-" << j + 1;
    } while (j != startNode);
    std::cout << std::endl;

    delete[] distance;
    delete[] visited;
    delete[] pred;
    delete[] pred2;
}



int main(int argc, char* argv[]) {
    const std::string fileWithSavedName = "N.txt";
    const std::string fileWithSavedPoint = "P.txt";
    const std::string fileWithDecision = "B.txt";
    const std::string fileWithEulerCycle = "E.txt";
    switches s;
    bool isFileLoaded = false;

    if (argc == 3) {
        for (int i = 1; i < argc; i++) {
            auto x = std::string(argv[i]);

            if (x == "-i") {
                s.i = argv[i + 1];
                auto v = LoadFromFile(s.i, isFileLoaded, fileWithDecision);
                std::string stringIsFileLoaded = LoadNew(fileWithDecision);

                if (stringIsFileLoaded == "1") {
                    SaveNew(s.i, fileWithSavedName);
                }
                PrintInputFile(v);
            }
            else if (x == "-o") {
                s.o = argv[i + 1];
                std::string outFileName = s.o;
                std::string stringIsFileLoaded = LoadNew(fileWithDecision);

                if (stringIsFileLoaded == "1") {
                    std::string currentInFileName = LoadNew(fileWithSavedName);
                    std::string stringStartPoint = LoadNew(fileWithSavedPoint);
                    int startPoint = std::stoi(stringStartPoint);
                    auto v = LoadFromFile(currentInFileName, isFileLoaded, fileWithDecision);
                    size_t vectorSize = v.size();
                    bool isStartPoint = IsStartPoint(startPoint, v);

                    if (isStartPoint) {
                        auto g = LoadGraph(v);
                        size_t graphSize = g.size();
                        auto matrix = LoadMatrix(g);

                        bool isGraphConnected = IsGraphConnected(vectorSize, graphSize, matrix);

                        if (isGraphConnected) {
                            int nrOfOddNodes = CountOddNodes(g);
                            std::string stringStartPoint = LoadNew(fileWithSavedPoint);
                            int startPoint = std::stoi(stringStartPoint);
                            startPoint--;

                            if (nrOfOddNodes == 0) {
                                std::ofstream out{ fileWithEulerCycle };
                                if (out) {
                                    FleuryAlgorithm(startPoint, vectorSize, graphSize, matrix, out);
                                }
                                out.close();

                                auto e = LoadEulerCycle(fileWithEulerCycle);
                                int iterator = 0;

                                std::ofstream fout{ outFileName };
                                if (fout) {
                                    std::cout << std::endl;
                                    MakeOutputFile(e, v, iterator, fout);
                                    std::cout << std::endl;
                                }
                                fout.close();
                            }
                            else if (nrOfOddNodes >= 2) {
                                auto matrix2 = LoadMatrix2(v, g);
                                auto n = LoadOddNodes(g, nrOfOddNodes);
                                auto l = CombinePairs(n);

                                for (int i = 0; i < l.size(); i++) {
                                     DijkstraAlgorithm(matrix2, l[i].first-1, l[i].second - 1);
                                 }

                                //int nrOfPairings = CountNrOfPairings(nrOfOddNodes);;
                                 //Find the pairing with minimum shortest path connecting pairs.
                                //Dla k  wierzchołków liczba skojarzeń wierzchołków w pary jest równa: S  = ( k  - 1 ) • ( k  - 3 ) • ( k  - 5 ) • ... • 7 • 5 • 3 • 1
                            }
                        }
                        else {
                            std::cout << std::endl;
                            std::cout << "Graf nie jest spojny!" << std::endl;
                            std::cout << std::endl;
                        }
                    }
                    else {
                        std::cout << std::endl;
                        std::cout << "Wybrany punkt startowy nie istnieje dla wybranego pliku wejsciowego." << std::endl;
                        std::cout << std::endl;
                    }
                }
                else {
                    std::cout << std::endl;
                    std::cout << "Nie wybrano pliku wejsciowego!" << std::endl;
                    std::cout << std::endl;
                }
            }
            else if (x == "-p") {
                s.p = std::atoi(argv[i + 1]);
                std::string stringIsFileLoaded = LoadNew(fileWithDecision);

                if (stringIsFileLoaded == "1") {
                    std::string currentInFileName = LoadNew(fileWithSavedName);
                    auto v = LoadFromFile(currentInFileName, isFileLoaded, fileWithDecision);
                    bool isStartPoint = IsStartPoint(s.p, v);
                    std::string stringStartPoint = std::to_string(s.p);
                    SaveNew(stringStartPoint, fileWithSavedPoint);

                    if (isStartPoint) {
                        std::cout << std::endl;
                        std::cout << "Wybrano punkt startowy - " << s.p << std::endl;
                        std::cout << std::endl;
                    }
                    else {
                        std::cout << std::endl;
                        std::cout << "Wybrany punkt startowy nie istnieje dla wybranego pliku wejsciowego." << std::endl;
                        std::cout << std::endl;
                    }
                }
                else {
                    std::cout << std::endl;
                    std::cout << "Nie wybrano pliku wejsciowego!" << std::endl;
                    std::cout << std::endl;
                }
            }
        }
    }
    else {
        std::cout << std::endl;
        std::cout << "Niepoprawna ilosc argumentow! Sproboj ponownie." << std::endl;
        std::cout << std::endl;
    }
}