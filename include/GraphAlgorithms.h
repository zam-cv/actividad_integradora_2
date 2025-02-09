/**
 * @file GraphAlgorithms.h
 * @brief Implementación de algoritmos de grafos
 * @version 1.0
 * @date 2024-02-07
 */
#pragma once

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <fstream>
#include <sstream>
#include <optional>
#include "Edge.h"
#include "DisjointSet.h"
#include "Coordinate.h"

/**
 * @struct GraphData
 * @brief Estructura que almacena todos los datos necesarios para los algoritmos de grafos
 */
struct GraphData {
    /**
     * @brief Número total de vértices en el grafo
     * @details Representa la cantidad de nodos o colonias en el sistema
     */
    int numVertices;

    /**
     * @brief Vector que contiene todas las aristas del grafo
     * @details Cada arista contiene dos vértices (node1, node2) y su peso correspondiente
     */
    std::vector<Edge> edges;

    /**
     * @brief Matriz de adyacencia del grafo
     * @details Matriz de tamaño numVertices x numVertices donde graph[i][j] representa
     * el peso de la arista entre los vértices i y j. Un valor de 0 indica que no hay conexión
     */
    std::vector<std::vector<int>> graph;

    /**
     * @brief Matriz de capacidades para el algoritmo de flujo máximo
     * @details Matriz que almacena las capacidades máximas de transmisión entre cada par de vértices
     * Utilizada en el algoritmo de Ford-Fulkerson
     */
    std::vector<std::vector<int>> capacityMatrix;

    /**
     * @brief Vector de coordenadas de las centrales eléctricas
     * @details Almacena las ubicaciones (x,y) de todas las centrales eléctricas en el sistema
     */
    std::vector<Coordinate> centrales;

    /**
     * @brief Coordenadas de la nueva casa a conectar
     * @details Almacena la ubicación (x,y) de la nueva casa que necesita ser conectada
     * a la central más cercana
     */
    Coordinate newHouse;

    /**
     * @brief Mapeo de índices numéricos a identificadores de colonias
     * @details Asocia cada índice de vértice (0,1,2,...) con una letra (A,B,C,...)
     * para identificar las colonias de manera más legible
     */
    std::map<int, char> colonyMap;
};

/*
 * Implementación del algoritmo de Kruskal
 * @param edges: Vector que contiene todas las aristas del grafo
 * @param numVertices: Número de vértices
 * @param colonyMap: Mapeo de índices a nombres de colonias
 * @return: Vector de aristas en el MST
 */
std::vector<Edge> kruskal_mst(std::vector<Edge>& edges, int numVertices, std::map<int, char>& colonyMap) {
    std::sort(edges.begin(), edges.end()); 
    DisjointSet ds(numVertices);
    std::vector<Edge> mstEdges;

    for (const Edge& edge : edges) {
        if (ds.FindSet(edge.node1) != ds.FindSet(edge.node2)) {
            ds.UnionSets(edge.node1, edge.node2);
            mstEdges.push_back(edge);
        }
    }

    std::cout << "\n1." << std::endl;
    for (const Edge& edge : mstEdges) {
        std::cout << "(" << colonyMap[edge.node1] << ", " << colonyMap[edge.node2] << ")" << std::endl;
    }

    return mstEdges;
}

/*
 * Implementación del TSP usando backtracking
 * @param graph: Matriz de adyacencia del grafo
 * @param path: Camino actual siendo explorado
 * @param visited: Vector para rastrear nodos visitados
 * @param pos: Posición actual en el grafo
 * @param count: Número de nodos visitados
 * @param cost: Costo actual del camino
 * @param minCost: Costo mínimo
 * @param bestPath: Mejor camino
 */
int tsp_backtracking(std::vector<std::vector<int>>& graph, std::vector<int>& path, 
                   std::vector<bool>& visited, int pos, int count, int cost, 
                   int& minCost, std::vector<int>& bestPath) {
    int n = graph.size();

    if (count == n && graph[pos][0] > 0) {
        int totalCost = cost + graph[pos][0];
        if (totalCost < minCost) {
            minCost = totalCost;
            bestPath = path;
            bestPath.push_back(0);
        }
        return minCost;
    }

    for (int i = 0; i < n; i++) {
        if (!visited[i] && graph[pos][i] > 0) {
            visited[i] = true;
            path.push_back(i);
            tsp_backtracking(graph, path, visited, i, count + 1, cost + graph[pos][i], minCost, bestPath);
            path.pop_back();
            visited[i] = false;
        }
    }

    return minCost;
}

/*
 * Resuelve el problema del TSP
 * @param graph: Matriz de adyacencia del grafo
 * @param colonyMap: Mapeo de índices a nombres de colonias
 */
void solve_tsp(std::vector<std::vector<int>>& graph, std::map<int, char>& colonyMap) {
    int n = graph.size();
    std::vector<int> path = {0}, bestPath;
    std::vector<bool> visited(n, false);
    visited[0] = true;
    int minCost = std::numeric_limits<int>::max();
    tsp_backtracking(graph, path, visited, 0, 1, 0, minCost, bestPath);

    std::cout << "\n2." << std::endl;
    for (int idx : bestPath) {
        std::cout << colonyMap[idx] << " ";
    }
    std::cout << std::endl;
}

/*
 * Realiza una búsqueda en profundidad para encontrar un camino de aumento
 * @param residualGraph: Grafo residual que muestra la capacidad restante
 * @param source: Nodo fuente
 * @param sink: Nodo sumidero
 * @param parent: Arreglo para almacenar el camino
 * @param visited: Arreglo para rastrear nodos visitados
 * @return: Verdadero si existe un camino, falso en caso contrario
 */
bool find_path(std::vector<std::vector<int>>& residualGraph, int source, int sink, 
             std::vector<int>& parent, std::vector<bool>& visited) {
    int numVertices = residualGraph.size();
    visited[source] = true;

    // Try all vertices as next step from source
    for (int v = 0; v < numVertices; v++) {
        // If not visited and has capacity
        if (!visited[v] && residualGraph[source][v] > 0) {
            parent[v] = source;

            // If we reached sink, path found
            if (v == sink) {
                return true;
            }

            // Recursively explore from v
            if (find_path(residualGraph, v, sink, parent, visited)) {
                return true;
            }
        }
    }
    return false;
}

/*
 * Implementación del algoritmo de Ford-Fulkerson
 * @param graph: Grafo original de capacidades
 * @param source: Nodo fuente
 * @param sink: Nodo sumidero
 * @return: Valor del flujo máximo
 */
int ford_fulkerson(std::vector<std::vector<int>>& graph, int source, int sink) {
    int numVertices = graph.size();
    std::vector<std::vector<int>> residualGraph = graph; // Create residual graph
    std::vector<int> parent(numVertices);
    int maxFlow = 0;

    // While there exists an augmenting path
    while (true) {
        std::vector<bool> visited(numVertices, false);

        // If no path exists, we're done
        if (!find_path(residualGraph, source, sink, parent, visited)) {
            break;
        }

        // Update residual capacities
        int pathFlow = std::numeric_limits<int>::max();
        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            pathFlow = std::min(pathFlow, residualGraph[u][v]);
        }

        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            residualGraph[u][v] -= pathFlow;
            residualGraph[v][u] += pathFlow;
        }

        maxFlow += pathFlow;
    }

    return maxFlow;
}

/*
 * Lee la matriz de capacidad y calcula el flujo máximo
 * @param capacityMatrix: Matriz con capacidades máximas de transmisión
 * @param colonyMap: Mapeo de índices a nombres de colonias
 */
void calculate_max_flow(std::vector<std::vector<int>>& capacityMatrix, std::map<int, char>& colonyMap) {
    int source = 0; // First colony
    int sink = capacityMatrix.size() - 1; // Last colony
    int maxFlow = ford_fulkerson(capacityMatrix, source, sink);
    std::cout << "\n3." << std::endl;
    std::cout << maxFlow << std::endl;
}

/*
 * Lee coordenadas desde un string en el formato especificado
 * @param input: String con las coordenadas en el formato especificado
 * @param newHouse: Referencia para almacenar las coordenadas de la nueva casa
 * @return: Vector de coordenadas
 */
std::vector<Coordinate> read_coordinates(const std::string& input, Coordinate& newHouse) {
    std::vector<Coordinate> coordinates;
    std::istringstream iss(input);
    std::string line;

    std::getline(iss, line);

    while (std::getline(iss, line)) {
        if (line.empty()) {
            continue;
        }

        line = line.substr(1, line.length() - 2);
        size_t commaPos = line.find(',');

        if (commaPos != std::string::npos) {
            Coordinate coord;
            coord.x = std::stoi(line.substr(0, commaPos));
            coord.y = std::stoi(line.substr(commaPos + 1));
            coordinates.push_back(coord);
        }
    }

    if (!coordinates.empty()) {
        newHouse = coordinates.back();
        coordinates.pop_back();  
    } else {
        std::cout << "Error: not coordinates found" << std::endl;
    }

    return coordinates;
}

/*
 * Encuentra la central más cercana a la nueva casa
 * @param coordinates: Lista de coordenadas de las centrales
 * @param newHouse: Coordenadas de la nueva casa
 * @return: Coordenadas de la central más cercana
 */
Coordinate find_nearest_central(const std::vector<Coordinate>& coordinates, Coordinate newHouse) {
    if (coordinates.empty()) {
        std::cout << "Error: not centrals found" << std::endl;
        return {0, 0}; // Return invalid coordinates
    }

    Coordinate nearestCentral;
    double minDistance = std::numeric_limits<double>::max();

    // Iterate through all centrals and compute Euclidean distance
    for (const auto& central : coordinates) {
        double distance = std::sqrt(std::pow(central.x - newHouse.x, 2) + 
                                  std::pow(central.y - newHouse.y, 2));
        if (distance < minDistance) {
            minDistance = distance;
            nearestCentral = central;
        }
    }

    return nearestCentral;
}

/*
 * Lee las aristas desde un archivo de entrada
 * @param inputFile: Flujo de entrada
 * @param numVertices: Número de vértices
 * @param graph: Matriz de adyacencia
 * @return: Vector de aristas
 */
std::vector<Edge> read_edges(std::istream& inputFile, int numVertices, std::vector<std::vector<int>>& graph) {
    std::vector<Edge> edges;
    for (int i = 0; i < numVertices; i++) {
        for (int j = 0; j < numVertices; j++) {
            inputFile >> graph[i][j];
            if (i < j && graph[i][j] > 0) {
                edges.push_back({i, j, graph[i][j]});
            }
        }
    }
    auto result = edges;
    return result;
}

/*
 * Lee la matriz de capacidad desde un archivo de entrada
 * @param inputFile: Flujo de entrada
 * @param numVertices: Número de vértices
 * @return: Matriz de capacidad
 */
std::vector<std::vector<int>> read_capacity_matrix(std::istream& inputFile, int numVertices) {
    std::vector<std::vector<int>> capacityMatrix(numVertices, std::vector<int>(numVertices));
    for (int i = 0; i < numVertices; i++) {
        for (int j = 0; j < numVertices; j++) {
            inputFile >> capacityMatrix[i][j];
        }
    }
    auto result = capacityMatrix;
    return result;
}

/*
 * Lee las coordenadas desde un archivo de entrada
 * @param inputFile: Flujo de entrada
 * @return: Cadena de coordenadas
 */
std::string read_coordinates_input(std::istream& inputFile) {
    std::string coordinatesInput;
    std::string line;
    std::getline(inputFile, line); // consume empty line
    while (std::getline(inputFile, line)) {
        coordinatesInput += line + "\n";
    }
    return coordinatesInput;
}

/*
 * Procesa el archivo de entrada y retorna una estructura con todos los datos necesarios
 * @param filename: Nombre del archivo de entrada
 * @return: tuple con todos los datos procesados, o nullopt si hay error
 */
std::optional<GraphData> process_input(std::istream& input) {
    GraphData data;
    input >> data.numVertices;
    
    data.graph.resize(data.numVertices, std::vector<int>(data.numVertices));
    for (int i = 0; i < data.numVertices; i++) {
        data.colonyMap[i] = 'A' + i;
    }

    data.edges = read_edges(input, data.numVertices, data.graph);
    data.capacityMatrix = read_capacity_matrix(input, data.numVertices);
    std::string coordinatesInput = read_coordinates_input(input);

    data.centrales = read_coordinates(coordinatesInput, data.newHouse);
    if (data.centrales.empty()) {
        return std::nullopt;
    }

    return data;
}

/**
 * Ejecuta todos los algoritmos de grafos y muestra los resultados
 * @param graphData: Estructura con todos los datos del grafo
 * @return: true si la ejecución fue exitosa, false en caso contrario
 */
bool execute_algorithms(const GraphData* graphData) {
    if (!graphData) {
        return false;
    }

    std::vector<Edge> edges = graphData->edges;
    std::vector<std::vector<int>> graph = graphData->graph;
    std::vector<std::vector<int>> capacityMatrix = graphData->capacityMatrix;
    std::map<int, char> colonyMap = graphData->colonyMap;

    kruskal_mst(edges, graphData->numVertices, colonyMap);
    solve_tsp(graph, colonyMap);
    calculate_max_flow(capacityMatrix, colonyMap);
    Coordinate nearestCentral = find_nearest_central(graphData->centrales, graphData->newHouse);
    std::cout << "\n4." << std::endl << "(" << nearestCentral.x << ", " << nearestCentral.y << ")" << std::endl;

    return true;
}

/**
 * @brief Procesa el archivo de entrada y ejecuta los algoritmos
 * @param inputFile Archivo de entrada
 * @return int Código de error
 */
int process_file_and_execute(std::ifstream& inputFile) {
    if (!inputFile) {
        std::cerr << "Error: No se pudo abrir el archivo de entrada" << std::endl;
        return 1;
    }

    auto graphData = process_input(inputFile);
    if (!graphData) {
        std::cerr << "Error: Datos de entrada inválidos" << std::endl;
        return 1;
    }
    
    execute_algorithms(&(*graphData));

    return 0;
}
