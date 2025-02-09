#include <gtest/gtest.h>
#include <vector>
#include <sstream>
#include <map>
#include "Edge.h"
#include "GraphAlgorithms.h"
#include "Coordinate.h"

TEST(graph_algorithms_test, kruskal_mst_basic) {
    std::vector<Edge> edges = {
        {0, 1, 4},
        {0, 2, 3},
        {1, 2, 1},
        {1, 3, 2},
        {2, 3, 4}
    };
    int numVertices = 4;
    std::map<int, char> colonyMap = {{0, 'A'}, {1, 'B'}, {2, 'C'}, {3, 'D'}};

    auto result = kruskal_mst(edges, numVertices, colonyMap);

    EXPECT_EQ(result.size(), numVertices - 1);

    int totalWeight = 0;
    for (const auto& edge : result) {
        totalWeight += edge.weight;
    }
    EXPECT_EQ(totalWeight, 6);

    std::set<int> vertices;
    for (const auto& edge : result) {
        vertices.insert(edge.node1);
        vertices.insert(edge.node2);
    }
    EXPECT_EQ(vertices.size(), numVertices);

    std::vector<Edge> expectedEdges = {
        {1, 2, 1},
        {1, 3, 2},
        {0, 2, 3}
    };
    
    std::sort(result.begin(), result.end());
    std::sort(expectedEdges.begin(), expectedEdges.end());
    EXPECT_EQ(result, expectedEdges);
}

TEST(graph_algorithms_test, kruskal_mst_empty) {
    std::vector<Edge> edges;
    int numVertices = 0;
    std::map<int, char> colonyMap;

    auto result = kruskal_mst(edges, numVertices, colonyMap);
    EXPECT_TRUE(result.empty());
}

TEST(graph_algorithms_test, kruskal_mst_single_vertex) {
    std::vector<Edge> edges;
    int numVertices = 1;
    std::map<int, char> colonyMap = {{0, 'A'}};

    auto result = kruskal_mst(edges, numVertices, colonyMap);
    EXPECT_TRUE(result.empty());
}

TEST(graph_algorithms_test, kruskal_mst_disconnected) {
    std::vector<Edge> edges = {
        {0, 1, 1},
        {2, 3, 2}
    };
    int numVertices = 4;
    std::map<int, char> colonyMap = {{0, 'A'}, {1, 'B'}, {2, 'C'}, {3, 'D'}};

    auto result = kruskal_mst(edges, numVertices, colonyMap);
    EXPECT_EQ(result.size(), 2);
    EXPECT_EQ(result[0].weight + result[1].weight, 3);
}

TEST(graph_algorithms_test, ford_fulkerson_no_path) {
    std::vector<std::vector<int>> graph = {
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}
    };
    
    int source = 0;
    int sink = 3;
    int maxFlow = ford_fulkerson(graph, source, sink);
    
    EXPECT_EQ(maxFlow, 0);
}

TEST(graph_algorithms_test, ford_fulkerson_single_path) {
    std::vector<std::vector<int>> graph = {
        {0, 5, 0, 0},
        {0, 0, 4, 0},
        {0, 0, 0, 3},
        {0, 0, 0, 0}
    };
    
    int source = 0;
    int sink = 3;
    int maxFlow = ford_fulkerson(graph, source, sink);
    
    EXPECT_EQ(maxFlow, 3);
}

TEST(graph_algorithms_test, ford_fulkerson_multiple_paths) {
    std::vector<std::vector<int>> graph = {
        {0, 10, 10, 0},
        {0, 0, 2, 10},
        {0, 0, 0, 10},
        {0, 0, 0, 0}
    };
    
    int source = 0;
    int sink = 3;
    int maxFlow = ford_fulkerson(graph, source, sink);
    
    EXPECT_EQ(maxFlow, 20);
    
    for (size_t i = 0; i < graph.size(); i++) {
        for (size_t j = 0; j < graph[i].size(); j++) {
            if (i != j) {
                EXPECT_GE(graph[i][j], 0);
                EXPECT_GE(graph[j][i], 0);
            }
        }
    }
}

TEST(graph_algorithms_test, find_nearest_central) {
    std::vector<Coordinate> centrals = {
        {0, 0},
        {3, 4},
        {-2, 1}
    };
    Coordinate newHouse = {1, 2};
    
    Coordinate nearest = find_nearest_central(centrals, newHouse);
    
    bool found = false;
    for (const auto& central : centrals) {
        if (central.x == nearest.x && central.y == nearest.y) {
            found = true;
            break;
        }
    }
    EXPECT_TRUE(found);
    
    double distanceToNearest = std::sqrt(std::pow(nearest.x - newHouse.x, 2) + 
                                       std::pow(nearest.y - newHouse.y, 2));
    
    for (const auto& central : centrals) {
        double distance = std::sqrt(std::pow(central.x - newHouse.x, 2) + 
                                  std::pow(central.y - newHouse.y, 2));
        
        EXPECT_GE(distance, distanceToNearest);
        
        if (std::abs(distance - distanceToNearest) < 1e-10) {
            EXPECT_TRUE(central.x == nearest.x && central.y == nearest.y);
        }
    }
}

TEST(graph_algorithms_test, find_nearest_central_empty) {
    std::vector<Coordinate> centrals;
    Coordinate newHouse = {1, 2};
    
    testing::internal::CaptureStdout();
    Coordinate result = find_nearest_central(centrals, newHouse);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_EQ(result.x, 0);
    EXPECT_EQ(result.y, 0);
    EXPECT_TRUE(output.find("Error: not centrals found") != std::string::npos);
}

TEST(graph_algorithms_test, find_nearest_central_same_distance) {
    std::vector<Coordinate> centrals = {
        {1, 1},
        {1, -1},
        {-1, 1},
        {-1, -1}
    };
    Coordinate newHouse = {0, 0};
    
    Coordinate nearest = find_nearest_central(centrals, newHouse);
    
    double expectedDistance = std::sqrt(2);
    double actualDistance = std::sqrt(std::pow(nearest.x - newHouse.x, 2) + 
                                    std::pow(nearest.y - newHouse.y, 2));
    
    EXPECT_DOUBLE_EQ(actualDistance, expectedDistance);
}

TEST(graph_algorithms_test, find_path_basic) {
    std::vector<std::vector<int>> residualGraph = {
        {0, 1, 1, 0},
        {0, 0, 1, 1},
        {0, 0, 0, 1},
        {0, 0, 0, 0}
    };
    std::vector<int> parent(4, -1);
    std::vector<bool> visited(4, false);
    
    bool hasPath = find_path(residualGraph, 0, 3, parent, visited);
    
    EXPECT_TRUE(hasPath);
    
    int current = 3;
    std::vector<int> path;
    while (current != -1) {
        path.push_back(current);
        if (current == 0) {
            break;
        }
        current = parent[current];
    }
    
    if (!path.empty()) {
        std::reverse(path.begin(), path.end());
        EXPECT_EQ(path.front(), 0);
        EXPECT_EQ(path.back(), 3);
        
        for (size_t i = 0; i < path.size() - 1; i++) {
            EXPECT_GT(residualGraph[path[i]][path[i + 1]], 0);
            EXPECT_TRUE(visited[path[i]]);
        }
    }
}

TEST(graph_algorithms_test, find_path_no_path) {
    std::vector<std::vector<int>> residualGraph = {
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}
    };
    std::vector<int> parent(4);
    std::vector<bool> visited(4, false);
    
    bool hasPath = find_path(residualGraph, 0, 3, parent, visited);
    
    EXPECT_FALSE(hasPath);
}

TEST(graph_algorithms_test, kruskal_mst_duplicate_edges) {
    std::vector<Edge> edges = {
        {0, 1, 4},
        {0, 1, 2},
        {1, 2, 3},
        {1, 2, 1}
    };
    int numVertices = 3;
    std::map<int, char> colonyMap = {{0, 'A'}, {1, 'B'}, {2, 'C'}};

    auto result = kruskal_mst(edges, numVertices, colonyMap);

    EXPECT_EQ(result.size(), numVertices - 1);
    
    std::vector<Edge> expectedEdges = {
        {1, 2, 1},
        {0, 1, 2}
    };
    
    std::sort(result.begin(), result.end());
    std::sort(expectedEdges.begin(), expectedEdges.end());
    EXPECT_EQ(result, expectedEdges);
}

TEST(graph_algorithms_test, tsp_backtracking_basic) {
    std::vector<std::vector<int>> graph = {
        {0, 10, 15, 20},
        {10, 0, 35, 25},
        {15, 35, 0, 30},
        {20, 25, 30, 0}
    };
    std::vector<int> path = {0};
    std::vector<bool> visited(4, false);
    visited[0] = true;
    std::vector<int> bestPath;
    int minCost = std::numeric_limits<int>::max();
    
    int result = tsp_backtracking(graph, path, visited, 0, 1, 0, minCost, bestPath);
    
    EXPECT_EQ(result, 80);
    EXPECT_EQ(bestPath.size(), 5);
    EXPECT_EQ(bestPath[0], 0);
    EXPECT_EQ(bestPath.back(), 0);
}

TEST(graph_algorithms_test, tsp_backtracking_small_graph) {
    std::vector<std::vector<int>> graph = {
        {0, 1, 2},
        {1, 0, 3},
        {2, 3, 0}
    };
    std::vector<int> path = {0};
    std::vector<bool> visited(3, false);
    visited[0] = true;
    std::vector<int> bestPath;
    int minCost = std::numeric_limits<int>::max();
    
    int result = tsp_backtracking(graph, path, visited, 0, 1, 0, minCost, bestPath);
    
    EXPECT_EQ(result, 6);
    EXPECT_EQ(bestPath.size(), 4);
}

TEST(graph_algorithms_test, tsp_backtracking_disconnected) {
    std::vector<std::vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 0},
        {0, 0, 0}
    };
    std::vector<int> path = {0};
    std::vector<bool> visited(3, false);
    visited[0] = true;
    std::vector<int> bestPath;
    int minCost = std::numeric_limits<int>::max();
    
    int result = tsp_backtracking(graph, path, visited, 0, 1, 0, minCost, bestPath);
    
    EXPECT_EQ(result, std::numeric_limits<int>::max());
    EXPECT_TRUE(bestPath.empty());
}

TEST(graph_algorithms_test, solve_tsp_basic) {
    std::vector<std::vector<int>> graph = {
        {0, 10, 15, 20},
        {10, 0, 35, 25},
        {15, 35, 0, 30},
        {20, 25, 30, 0}
    };
    std::map<int, char> colonyMap = {{0, 'A'}, {1, 'B'}, {2, 'C'}, {3, 'D'}};
    
    testing::internal::CaptureStdout();
    solve_tsp(graph, colonyMap);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_TRUE(output.find("2.") != std::string::npos);
    EXPECT_TRUE(output.find("A") != std::string::npos);
    size_t firstA = output.find("A");
    size_t lastA = output.rfind("A");
    EXPECT_TRUE(firstA != std::string::npos);
    EXPECT_TRUE(lastA != std::string::npos);
    EXPECT_TRUE(firstA != lastA);
}

TEST(graph_algorithms_test, solve_tsp_small_circuit) {
    std::vector<std::vector<int>> graph = {
        {0, 1, 1},
        {1, 0, 1},
        {1, 1, 0}
    };
    std::map<int, char> colonyMap = {{0, 'A'}, {1, 'B'}, {2, 'C'}};
    
    testing::internal::CaptureStdout();
    solve_tsp(graph, colonyMap);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_TRUE(output.find("2.") != std::string::npos);
    EXPECT_TRUE(output.find("A") != std::string::npos);
    EXPECT_TRUE(output.find("B") != std::string::npos);
    EXPECT_TRUE(output.find("C") != std::string::npos);
}

TEST(graph_algorithms_test, solve_tsp_disconnected) {
    std::vector<std::vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 0},
        {0, 0, 0}
    };
    std::map<int, char> colonyMap = {{0, 'A'}, {1, 'B'}, {2, 'C'}};
    
    testing::internal::CaptureStdout();
    solve_tsp(graph, colonyMap);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_TRUE(output.find("2.") != std::string::npos);
}

TEST(graph_algorithms_test, solve_tsp_edge_cases) {
    std::vector<std::vector<int>> graph = {
        {0, std::numeric_limits<int>::max(), 1},
        {std::numeric_limits<int>::max(), 0, 1},
        {1, 1, 0}
    };
    std::map<int, char> colonyMap = {{0, 'A'}, {1, 'B'}, {2, 'C'}};
    
    testing::internal::CaptureStdout();
    solve_tsp(graph, colonyMap);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_TRUE(output.find("2.") != std::string::npos);

    for (const auto& pair : colonyMap) {
        EXPECT_TRUE(output.find(pair.second) != std::string::npos);
    }
}

TEST(graph_algorithms_test, calculate_max_flow_basic) {
    std::vector<std::vector<int>> capacityMatrix = {
        {0, 10, 10, 0},
        {0, 0, 2, 10},
        {0, 0, 0, 10},
        {0, 0, 0, 0}
    };
    std::map<int, char> colonyMap = {{0, 'A'}, {1, 'B'}, {2, 'C'}, {3, 'D'}};
    
    testing::internal::CaptureStdout();
    calculate_max_flow(capacityMatrix, colonyMap);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_TRUE(output.find("3.") != std::string::npos);
    EXPECT_TRUE(output.find("20") != std::string::npos);
}

TEST(graph_algorithms_test, calculate_max_flow_no_path) {
    std::vector<std::vector<int>> capacityMatrix = {
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}
    };
    std::map<int, char> colonyMap = {{0, 'A'}, {1, 'B'}, {2, 'C'}, {3, 'D'}};
    
    testing::internal::CaptureStdout();
    calculate_max_flow(capacityMatrix, colonyMap);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_TRUE(output.find("3.") != std::string::npos);
    EXPECT_TRUE(output.find("0") != std::string::npos);
}

TEST(graph_algorithms_test, calculate_max_flow_single_path) {
    std::vector<std::vector<int>> capacityMatrix = {
        {0, 5, 0, 0},
        {0, 0, 4, 0},
        {0, 0, 0, 3},
        {0, 0, 0, 0}
    };
    std::map<int, char> colonyMap = {{0, 'A'}, {1, 'B'}, {2, 'C'}, {3, 'D'}};
    
    testing::internal::CaptureStdout();
    calculate_max_flow(capacityMatrix, colonyMap);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_TRUE(output.find("3.") != std::string::npos);
    EXPECT_TRUE(output.find("3") != std::string::npos);
}

TEST(graph_algorithms_test, read_coordinates_basic) {
    std::string input = 
        "4\n"
        "(1,2)\n"
        "(3,4)\n"
        "(5,6)\n"
        "(7,8)\n";
    
    Coordinate newHouse;
    auto coordinates = read_coordinates(input, newHouse);
    
    EXPECT_EQ(coordinates.size(), 3);
    EXPECT_EQ(coordinates[0].x, 1);
    EXPECT_EQ(coordinates[0].y, 2);
    EXPECT_EQ(coordinates[1].x, 3);
    EXPECT_EQ(coordinates[1].y, 4);
    EXPECT_EQ(coordinates[2].x, 5);
    EXPECT_EQ(coordinates[2].y, 6);
    EXPECT_EQ(newHouse.x, 7);
    EXPECT_EQ(newHouse.y, 8);
}

TEST(graph_algorithms_test, read_coordinates_empty) {
    std::string input = "0\n";
    
    Coordinate newHouse;
    testing::internal::CaptureStdout();
    auto coordinates = read_coordinates(input, newHouse);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_TRUE(coordinates.empty());
    EXPECT_TRUE(output.find("Error: not coordinates found") != std::string::npos);
}

TEST(graph_algorithms_test, read_coordinates_single) {
    std::string input = 
        "1\n"
        "(1,2)\n";
    
    Coordinate newHouse;
    auto coordinates = read_coordinates(input, newHouse);
    
    EXPECT_TRUE(coordinates.empty());
    EXPECT_EQ(newHouse.x, 1);
    EXPECT_EQ(newHouse.y, 2);
}

TEST(graph_algorithms_test, read_coordinates_invalid_format) {
    std::string input = 
        "3\n"
        "(1,2)\n"
        "invalid_format\n"
        "(5,6)\n";
    
    Coordinate newHouse;
    auto coordinates = read_coordinates(input, newHouse);
    
    EXPECT_EQ(coordinates.size(), 1);
    EXPECT_EQ(coordinates[0].x, 1);
    EXPECT_EQ(coordinates[0].y, 2);
}

TEST(graph_algorithms_test, read_coordinates_empty_lines) {
    std::string input = 
        "4\n"
        "(1,2)\n"
        "\n"
        "(3,4)\n"
        "(5,6)\n";
    
    Coordinate newHouse;
    auto coordinates = read_coordinates(input, newHouse);
    
    EXPECT_EQ(coordinates.size(), 2);
    EXPECT_EQ(newHouse.x, 5);
    EXPECT_EQ(newHouse.y, 6);
}

TEST(graph_algorithms_test, read_coordinates_missing_comma) {
    std::string input = 
        "3\n"
        "(12)\n"
        "(3,4)\n"
        "(5,6)\n";
    
    Coordinate newHouse;
    auto coordinates = read_coordinates(input, newHouse);
    
    EXPECT_EQ(coordinates.size(), 1);
    EXPECT_EQ(coordinates[0].x, 3);
    EXPECT_EQ(coordinates[0].y, 4);
}

TEST(input_processing_test, read_edges_basic) {
    std::stringstream ss;
    ss << "0 1 2\n"
       << "1 0 3\n"
       << "2 3 0\n";
    std::vector<std::vector<int>> graph(3, std::vector<int>(3));
    
    auto edges = read_edges(ss, 3, graph);
    
    EXPECT_EQ(edges.size(), 3);
    EXPECT_EQ(edges[0].node1, 0);
    EXPECT_EQ(edges[0].node2, 1);
    EXPECT_EQ(edges[0].weight, 1);
}

TEST(input_processing_test, read_capacity_matrix_basic) {
    std::stringstream ss;
    ss << "0 5 3\n"
       << "0 0 2\n"
       << "0 0 0\n";
    
    auto matrix = read_capacity_matrix(ss, 3);
    
    EXPECT_EQ(matrix.size(), 3);
    EXPECT_EQ(matrix[0][1], 5);
    EXPECT_EQ(matrix[1][2], 2);
    EXPECT_EQ(matrix[0][0], 0);
}

TEST(input_processing_test, read_edges_empty) {
    std::stringstream ss;
    ss << "0 0 0\n"
       << "0 0 0\n"
       << "0 0 0\n";
    std::vector<std::vector<int>> graph(3, std::vector<int>(3));
    
    auto edges = read_edges(ss, 3, graph);
    
    EXPECT_TRUE(edges.empty());
}

TEST(input_processing_test, read_capacity_matrix_zeros) {
    std::stringstream ss;
    ss << "0 0 0\n"
       << "0 0 0\n"
       << "0 0 0\n";
    
    auto matrix = read_capacity_matrix(ss, 3);
    
    for (const auto& row : matrix) {
        for (int val : row) {
            EXPECT_EQ(val, 0);
        }
    }
}

TEST(input_processing_test, read_coordinates_input_empty) {
    std::stringstream ss;
    
    std::string result = read_coordinates_input(ss);
    
    EXPECT_TRUE(result.empty());
}

TEST(input_processing_test, process_input_basic) {
    std::stringstream ss;
    ss << "3\n"
       << "0 1 2\n"
       << "1 0 3\n"
       << "2 3 0\n"
       << "0 5 3\n"
       << "5 0 2\n"
       << "3 2 0\n"
       << "\n"
       << "(1,2)\n"
       << "(3,4)\n"
       << "(5,6)\n";

    auto result = process_input(ss);
    
    EXPECT_TRUE(result.has_value());
    EXPECT_EQ(result->numVertices, 3);
    EXPECT_EQ(result->edges.size(), 3);
    EXPECT_EQ(result->centrales.size(), 2);
    EXPECT_EQ(result->newHouse.x, 5);
    EXPECT_EQ(result->newHouse.y, 6);
    EXPECT_EQ(result->colonyMap[0], 'A');
    EXPECT_EQ(result->colonyMap[1], 'B');
    EXPECT_EQ(result->colonyMap[2], 'C');
}

TEST(input_processing_test, process_input_no_coordinates) {
    std::stringstream ss;
    ss << "3\n"
       << "0 1 2\n"
       << "1 0 3\n"
       << "2 3 0\n"
       << "0 5 3\n"
       << "5 0 2\n"
       << "3 2 0\n"
       << "\n";

    auto result = process_input(ss);
    
    EXPECT_FALSE(result.has_value());
}

TEST(input_processing_test, read_edges_negative_weights) {
    std::stringstream ss;
    ss << "0 -1 2\n"
       << "-1 0 -3\n"
       << "2 -3 0\n";
    std::vector<std::vector<int>> graph(3, std::vector<int>(3));
    
    auto edges = read_edges(ss, 3, graph);

    EXPECT_EQ(edges.size(), 1);
    EXPECT_EQ(edges[0].node1, 0);
    EXPECT_EQ(edges[0].node2, 2);
    EXPECT_EQ(edges[0].weight, 2);
}

TEST(input_processing_test, read_edges_large_graph) {
    std::stringstream ss;
    const int size = 5;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            ss << (i < j ? i + j : 0) << " ";
        }
        ss << "\n";
    }
    std::vector<std::vector<int>> graph(size, std::vector<int>(size));
    
    auto edges = read_edges(ss, size, graph);
    
    EXPECT_EQ(edges.size(), (size * (size - 1)) / 2);
    for (const auto& edge : edges) {
        EXPECT_LT(edge.node1, edge.node2);
        EXPECT_EQ(edge.weight, edge.node1 + edge.node2);
    }
}

TEST(input_processing_test, read_capacity_matrix_max_values) {
    std::stringstream ss;
    ss << "0 " << std::numeric_limits<int>::max() << " 3\n"
       << "0 0 " << std::numeric_limits<int>::max() << "\n"
       << "0 0 0\n";
    
    auto matrix = read_capacity_matrix(ss, 3);
    
    EXPECT_EQ(matrix[0][1], std::numeric_limits<int>::max());
    EXPECT_EQ(matrix[1][2], std::numeric_limits<int>::max());
}

TEST(input_processing_test, read_capacity_matrix_asymmetric) {
    std::stringstream ss;
    ss << "0 5 3\n"
       << "2 0 4\n"
       << "1 6 0\n";
    
    auto matrix = read_capacity_matrix(ss, 3);
    
    EXPECT_EQ(matrix[0][1], 5);
    EXPECT_EQ(matrix[1][0], 2);
    EXPECT_EQ(matrix[1][2], 4);
    EXPECT_EQ(matrix[2][1], 6);
}

TEST(input_processing_test, read_coordinates_input_multiple_empty_lines) {
    std::stringstream ss;
    ss << "\n\n(1,2)\n\n(3,4)\n\n";
    
    std::string result = read_coordinates_input(ss);
    
    EXPECT_TRUE(result.find("(1,2)") != std::string::npos);
    EXPECT_TRUE(result.find("(3,4)") != std::string::npos);
}

TEST(input_processing_test, process_input_large_data) {
    std::stringstream ss;
    const int size = 10;
    ss << size << "\n";
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            ss << (i < j ? 1 : 0) << " ";
        }
        ss << "\n";
    }
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            ss << (i < j ? 5 : 0) << " ";
        }
        ss << "\n";
    }
    
    ss << "\n";
    for (int i = 0; i < size; i++) {
        ss << "(" << i << "," << i << ")\n";
    }
    
    auto result = process_input(ss);
    
    EXPECT_TRUE(result.has_value());
    EXPECT_EQ(result->numVertices, size);
    EXPECT_EQ(result->edges.size(), (size * (size - 1)) / 2);
    EXPECT_EQ(result->centrales.size(), size - 1);
    EXPECT_EQ(result->colonyMap.size(), size);
}

TEST(input_processing_test, process_input_incomplete_data) {
    std::stringstream ss;
    ss << "3\n"
       << "0 1 2\n"
       << "1 0 3\n"
       << "2 3 0\n";
    
    auto result = process_input(ss);
    
    EXPECT_FALSE(result.has_value());
}

TEST(graph_algorithms_test, process_input_invalid_coordinates) {
    std::stringstream ss;
    ss << "3\n"
       << "0 1 2\n"
       << "1 0 3\n"
       << "2 3 0\n"
       << "0 5 3\n"
       << "5 0 2\n"
       << "3 2 0\n"
       << "\n"
       << "invalid_coordinate_format\n"
       << "(3,4)\n"
       << "(5,6)\n";

    auto result = process_input(ss);
    EXPECT_TRUE(result.has_value());
    EXPECT_EQ(result->centrales.size(), 1);
}

TEST(algorithm_execution_test, execute_algorithms_null_data) {
    EXPECT_FALSE(execute_algorithms(nullptr));
}

TEST(algorithm_execution_test, execute_algorithms_complex_data) {
    GraphData data;
    data.numVertices = 4;

    data.edges = {
        {0, 1, 16},
        {0, 2, 45},
        {0, 3, 32},
        {1, 2, 18},
        {1, 3, 21},
        {2, 3, 7}
    };
    
    data.graph = {
        {0, 16, 45, 32},
        {16, 0, 18, 21},
        {45, 18, 0, 7},
        {32, 21, 7, 0}
    };
    
    data.capacityMatrix = {
        {0, 48, 12, 18},
        {52, 0, 42, 32},
        {18, 46, 0, 56},
        {24, 36, 52, 0}
    };
    
    data.centrales = {
        {200, 500},
        {300, 100},
        {450, 150},
        {520, 480}
    };
    data.newHouse = {400, 300};
    
    data.colonyMap = {
        {0, 'A'},
        {1, 'B'},
        {2, 'C'},
        {3, 'D'}
    };

    testing::internal::CaptureStdout();
    EXPECT_TRUE(execute_algorithms(&data));
    std::string output = testing::internal::GetCapturedStdout();
    
    std::string expectedOutput = 
        "\n1.\n"
        "(C, D)\n"
        "(A, B)\n"
        "(B, C)\n"
        "\n"
        "2.\n"
        "A B C D A \n"
        "\n"
        "3.\n"
        "78\n"
        "\n"
        "4.\n"
        "(450, 150)\n";

    EXPECT_EQ(output, expectedOutput);
}

TEST(additional_tests, kruskal_mst_stdout_test) {
    std::vector<Edge> edges = { {0, 1, 4}, {0, 2, 3}, {1, 2, 1}, {1, 3, 2}, {2, 3, 4} };
    int numVertices = 4;
    std::map<int, char> colonyMap = { {0, 'A'}, {1, 'B'}, {2, 'C'}, {3, 'D'} };
    
    testing::internal::CaptureStdout();
    kruskal_mst(edges, numVertices, colonyMap);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_NE(output.find("\n1."), std::string::npos);
    EXPECT_NE(output.find("(A, "), std::string::npos);
}

TEST(additional_tests, solve_tsp_no_tour_test) {
    std::vector<std::vector<int>> graph = {
        {0, 10, 15},
        {10, 0, 20},
        {0, 20, 0}
    };
    std::map<int, char> colonyMap = { {0, 'A'}, {1, 'B'}, {2, 'C'} };
    
    testing::internal::CaptureStdout();
    solve_tsp(graph, colonyMap);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_NE(output.find("\n2."), std::string::npos);
}

TEST(additional_tests, calculate_max_flow_stdout_test) {
    std::vector<std::vector<int>> capacityMatrix = {
        {0, 10, 5, 0},
        {0, 0, 15, 10},
        {0, 0, 0, 10},
        {0, 0, 0, 0}
    };
    std::map<int, char> colonyMap = { {0, 'A'}, {1, 'B'}, {2, 'C'}, {3, 'D'} };
    
    testing::internal::CaptureStdout();
    calculate_max_flow(capacityMatrix, colonyMap);
    std::string output = testing::internal::GetCapturedStdout();
    
    EXPECT_NE(output.find("\n3."), std::string::npos);
}

TEST(additional_tests, read_coordinates_multiple_spaces_test) {
    std::string input = "3\n(1,2)\n   \n(3,4)\n(5,6)\n";
    Coordinate newHouse;
    auto coordinates = read_coordinates(input, newHouse);
    
    EXPECT_EQ(newHouse.x, 5);
    EXPECT_EQ(newHouse.y, 6);
    EXPECT_EQ(coordinates.size(), 2);
    EXPECT_EQ(coordinates[0].x, 1);
    EXPECT_EQ(coordinates[0].y, 2);
    EXPECT_EQ(coordinates[1].x, 3);
    EXPECT_EQ(coordinates[1].y, 4);
}

TEST(additional_tests, kruskal_no_valid_edges) {
    std::vector<Edge> edges = { {0, 1, 0}, {0, 2, 0} };
    int numVertices = 3;
    std::map<int, char> colonyMap = { {0, 'A'}, {1, 'B'}, {2, 'C'} };
    testing::internal::CaptureStdout();
    auto mst = kruskal_mst(edges, numVertices, colonyMap);
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(mst.size(), 2);
    EXPECT_NE(output.find("\n1."), std::string::npos);
}

TEST(additional_tests, tsp_no_available_moves) {
    std::vector<std::vector<int>> graph = {
        {0, 0},
        {0, 0}
    };
    std::vector<int> path = {0};
    std::vector<bool> visited(2, false);
    visited[0] = true;
    std::vector<int> bestPath;
    int minCost = std::numeric_limits<int>::max();
    int cost = tsp_backtracking(graph, path, visited, 0, 1, 0, minCost, bestPath);
    EXPECT_EQ(cost, std::numeric_limits<int>::max());
    EXPECT_TRUE(bestPath.empty());
}

TEST(additional_tests, find_path_source_equals_sink) {
    std::vector<std::vector<int>> residualGraph = {
        {0, 1},
        {0, 0}
    };
    std::vector<int> parent(2, -1);
    std::vector<bool> visited(2, false);

    bool result = find_path(residualGraph, 0, 0, parent, visited);
    EXPECT_FALSE(result);
}

TEST(additional_tests, read_edges_all_zeros) {
    std::stringstream ss;
    ss << "0 0\n0 0\n";
    int numVertices = 2;
    std::vector<std::vector<int>> graph(numVertices, std::vector<int>(numVertices, 0));
    auto edges = read_edges(ss, numVertices, graph);
    EXPECT_TRUE(edges.empty());
}

TEST(additional_tests, process_input_zero_vertices) {
    std::stringstream ss;
    ss << "0\n";
    auto result = process_input(ss);
    EXPECT_FALSE(result.has_value());
}

TEST(additional_tests, find_nearest_central_far_distance) {
    std::vector<Coordinate> centrals = { {100, 100}, {200, 200}, {300, 300} };
    Coordinate newHouse = {0, 0};
    testing::internal::CaptureStdout();
    Coordinate nearest = find_nearest_central(centrals, newHouse);
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(nearest.x, 100);
    EXPECT_EQ(nearest.y, 100);
}

TEST(additional_tests, read_capacity_matrix_negative_values) {
    std::stringstream ss;
    ss << "-1 -2\n-3 -4\n";
    auto matrix = read_capacity_matrix(ss, 2);
    EXPECT_EQ(matrix[0][0], -1);
    EXPECT_EQ(matrix[0][1], -2);
    EXPECT_EQ(matrix[1][0], -3);
    EXPECT_EQ(matrix[1][1], -4);
}

TEST(additional_tests, read_coordinates_extra_whitespace) {
    std::string input = "2\n( 7 , 8 )\n( 9 , 10 )\n";
    Coordinate newHouse;
    testing::internal::CaptureStdout();
    auto coordinates = read_coordinates(input, newHouse);
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(newHouse.x, 9);
    EXPECT_EQ(newHouse.y, 10);
    ASSERT_EQ(coordinates.size(), 1);
    EXPECT_EQ(coordinates[0].x, 7);
    EXPECT_EQ(coordinates[0].y, 8);
}

TEST(process_file_test, file_not_found) {
    std::ifstream nonexistentFile("nonexistent.txt");
    testing::internal::CaptureStderr();
    
    int result = process_file_and_execute(nonexistentFile);
    
    std::string error_output = testing::internal::GetCapturedStderr();
    EXPECT_EQ(result, 1);
    EXPECT_EQ(error_output, "Error: No se pudo abrir el archivo de entrada\n");
}

TEST(process_file_test, invalid_input_data) {
    std::string tempFileName = "temp_test.txt";
    std::ofstream tempFile(tempFileName);
    tempFile << "0\n";
    tempFile.close();
    
    std::ifstream inputFile(tempFileName);
    testing::internal::CaptureStderr();
    
    int result = process_file_and_execute(inputFile);
    
    std::string error_output = testing::internal::GetCapturedStderr();
    EXPECT_EQ(result, 1);
    EXPECT_EQ(error_output, "Error: Datos de entrada inválidos\n");
    
    inputFile.close();
    std::remove(tempFileName.c_str());
}

TEST(process_file_test, missing_coordinates) {
    std::string tempFileName = "temp_test.txt";
    std::ofstream tempFile(tempFileName);

    tempFile << "3\n"
             << "0 1 2\n"
             << "1 0 3\n"
             << "2 3 0\n"
             << "0 5 3\n"
             << "5 0 2\n"
             << "3 2 0\n"
             << "\n";
    tempFile.close();
    
    std::ifstream inputFile(tempFileName);
    testing::internal::CaptureStderr();
    
    int result = process_file_and_execute(inputFile);
    
    std::string error_output = testing::internal::GetCapturedStderr();
    EXPECT_EQ(result, 1);
    EXPECT_EQ(error_output, "Error: Datos de entrada inválidos\n");
    
    inputFile.close();
    std::remove(tempFileName.c_str());
}

TEST(process_file_test, valid_input) {
    std::string tempFileName = "temp_test.txt";
    std::ofstream tempFile(tempFileName);
    tempFile << "3\n"
             << "0 1 2\n"
             << "1 0 3\n"
             << "2 3 0\n"
             << "0 5 3\n"
             << "5 0 2\n"
             << "3 2 0\n"
             << "\n"
             << "(1,2)\n"
             << "(3,4)\n"
             << "(5,6)\n"
             << "(7,8)\n";
    tempFile.close();
    
    std::ifstream inputFile(tempFileName);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    
    int result = process_file_and_execute(inputFile);
    
    std::string stdout_output = testing::internal::GetCapturedStdout();
    std::string stderr_output = testing::internal::GetCapturedStderr();
    
    EXPECT_EQ(result, 0);
    EXPECT_TRUE(stderr_output.empty());
    EXPECT_FALSE(stdout_output.empty());

    EXPECT_NE(stdout_output.find("\n1."), std::string::npos);
    EXPECT_NE(stdout_output.find("\n2."), std::string::npos);
    EXPECT_NE(stdout_output.find("\n3."), std::string::npos);
    EXPECT_NE(stdout_output.find("\n4."), std::string::npos);
    
    inputFile.close();
    std::remove(tempFileName.c_str());
}

TEST(process_file_test, incomplete_input) {
    std::string tempFileName = "temp_test.txt";
    std::ofstream tempFile(tempFileName);
    tempFile << "3\n"
             << "0 1 2\n"
             << "1 0 3\n";
    tempFile.close();
    
    std::ifstream inputFile(tempFileName);
    testing::internal::CaptureStderr();
    
    int result = process_file_and_execute(inputFile);
    
    std::string error_output = testing::internal::GetCapturedStderr();
    EXPECT_EQ(result, 1);
    EXPECT_EQ(error_output, "Error: Datos de entrada inválidos\n");
    
    inputFile.close();
    std::remove(tempFileName.c_str());
}

TEST(process_file_test, invalid_matrix_format) {
    std::string tempFileName = "temp_test.txt";
    std::ofstream tempFile(tempFileName);

    tempFile << "2\n"
             << "0 a\n"
             << "b 0\n"
             << "0 1\n"
             << "1 0\n"
             << "\n"
             << "(1,2)\n"
             << "(3,4)\n"
             << "(5,6)\n";
    tempFile.close();
    
    std::ifstream inputFile(tempFileName);
    testing::internal::CaptureStderr();
    
    int result = process_file_and_execute(inputFile);
    
    std::string error_output = testing::internal::GetCapturedStderr();
    EXPECT_EQ(result, 1);
    EXPECT_EQ(error_output, "Error: Datos de entrada inválidos\n");
    
    inputFile.close();
    std::remove(tempFileName.c_str());
}

TEST(input_processing_test, read_edges_basic_with_return) {
    std::stringstream ss;
    ss << "0 1 2\n"
       << "1 0 3\n"
       << "2 3 0\n";
    std::vector<std::vector<int>> graph(3, std::vector<int>(3));
    
    auto edges = read_edges(ss, 3, graph);
    
    EXPECT_EQ(edges.size(), 3);
    
    bool foundEdge1 = false;
    bool foundEdge2 = false;
    bool foundEdge3 = false;
    
    for (const auto& edge : edges) {
        if (edge.node1 == 0 && edge.node2 == 1 && edge.weight == 1) foundEdge1 = true;
        if (edge.node1 == 1 && edge.node2 == 2 && edge.weight == 3) foundEdge2 = true;
        if (edge.node1 == 0 && edge.node2 == 2 && edge.weight == 2) foundEdge3 = true;
    }
    
    EXPECT_TRUE(foundEdge1);
    EXPECT_TRUE(foundEdge2);
    EXPECT_TRUE(foundEdge3);
}

TEST(input_processing_test, read_capacity_matrix_basic_with_return) {
    std::stringstream ss;
    ss << "0 5 3\n"
       << "0 0 2\n"
       << "0 0 0\n";
    
    auto matrix = read_capacity_matrix(ss, 3);

    std::vector<std::vector<int>> expected = {
        {0, 5, 3},
        {0, 0, 2},
        {0, 0, 0}
    };
    
    EXPECT_EQ(matrix.size(), expected.size());
    for (size_t i = 0; i < matrix.size(); i++) {
        EXPECT_EQ(matrix[i].size(), expected[i].size());
        for (size_t j = 0; j < matrix[i].size(); j++) {
            EXPECT_EQ(matrix[i][j], expected[i][j]);
        }
    }
}
