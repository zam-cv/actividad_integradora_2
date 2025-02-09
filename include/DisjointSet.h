/**
 * @file DisjointSet.h
 * @brief Implementación de la estructura de datos Disjoint Set
 * @version 1.0
 * @date 2024-02-07
 */

#pragma once
#include <vector>

/**
 * @class DisjointSet
 * @brief Clase para manejar conjuntos disjuntos con compresión de ruta y unión por rango
 */
class DisjointSet {
public:
    /** @brief Vector que almacena el padre de cada nodo */
    std::vector<int> parent;
    
    /** @brief Vector que almacena el rango (profundidad) de cada árbol */
    std::vector<int> rank;

    /**
     * @brief Constructor de la clase DisjointSet
     * @param n Número de elementos en el conjunto disjunto
     */
    DisjointSet(int n) {
        parent.resize(n);
        rank.resize(n, 0);
        for (int i = 0; i < n; i++) {
            parent[i] = i; 
        }
    }

    /**
     * @brief Encuentra el representante del conjunto al que pertenece un nodo
     * @param node Nodo a buscar
     * @return Representante del conjunto
     */
    int FindSet(int node) {
        if (parent[node] != node) {
            parent[node] = FindSet(parent[node]); //Path compression
        }
        return parent[node];
    }

    /**
     * @brief Une dos conjuntos disjuntos usando unión por rango
     * @param node1 Primer nodo a unir
     * @param node2 Segundo nodo a unir
     */
    void UnionSets(int node1, int node2) {
        int root1 = FindSet(node1);
        int root2 = FindSet(node2);
        if (root1 != root2) {
            if (rank[root1] < rank[root2]) {
                parent[root1] = root2;
            } else if (rank[root1] > rank[root2]) {
                parent[root2] = root1;
            } else {
                parent[root2] = root1;
                rank[root1]++;
            }
        }
    }
};
