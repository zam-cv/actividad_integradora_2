/**
 * @file Edge.h
 * @brief Define la estructura para representar aristas en un grafo
 * @version 1.0
 * @date 2024-02-07
 */

#pragma once

/**
 * @struct Edge
 * @brief Estructura para representar una arista en un grafo
 */
struct Edge {
    /** @brief Primer nodo de la arista */
    int node1;    
    /** @brief Segundo nodo de la arista */
    int node2;    
    /** @brief Peso/costo de la arista */
    int weight;   
    
    /**
     * @brief Operador de comparación para ordenar aristas por peso
     * @param other La arista con la que comparar
     * @return true si el peso de esta arista es menor que el de la otra
     */
    bool operator<(const Edge& other) const {
        if (weight != other.weight) {
            return weight < other.weight;
        }
        if (node1 != other.node1) {
            return node1 < other.node1;
        }
        return node2 < other.node2;
    }

    /**
     * @brief Operador de igualdad para comparar aristas
     * @param other La arista con la que comparar
     * @return true si las aristas son iguales
     */
    bool operator==(const Edge& other) const {
        return node1 == other.node1 && 
               node2 == other.node2 && 
               weight == other.weight;
    }
};
