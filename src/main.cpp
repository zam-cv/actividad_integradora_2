/*
 * Actividad Integradora 2
 * 
 * Program that implements various string pattern matching algorithms:
 * - (Kruskal's)
 * - (Traveling salesman)
 * - (Ford Fulkenson)
 * - (Búsqueda lineal)
 * 
 * Authors:
 * - Lorena Abigail Solís de los Santos   A01746602
 * - Andrea Doce Murillo                  A01799931
 * - Carlos Alberto Zamudio Velázquez     A01799283
 * 
 * Date: 07/02/2025
 */
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <map>
#include <limits>
#include "Edge.h"
#include "DisjointSet.h"
#include "Coordinate.h"
#include "GraphAlgorithms.h"

using namespace std;

int main() {
    std::ifstream inputFile("input.txt");
    return process_file_and_execute(inputFile);
}
