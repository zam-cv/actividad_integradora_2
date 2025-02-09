#include <gtest/gtest.h>
#include <climits>
#include "Edge.h"

class EdgeTest : public ::testing::Test {
protected:
    void SetUp() override {
        edge1 = Edge{1, 2, 10};
        edge2 = Edge{1, 2, 20};
        edge3 = Edge{1, 2, 10};
        edge4 = Edge{2, 3, 5};
    }

    Edge edge1, edge2, edge3, edge4;
};

TEST_F(EdgeTest, OperadorMenorQue) {
    EXPECT_TRUE(edge1 < edge2);
    EXPECT_FALSE(edge2 < edge1);
    EXPECT_FALSE(edge1 < edge3);

    Edge edge5{0, 2, 10};
    Edge edge6{3, 2, 10};
    EXPECT_TRUE(edge5 < edge1);
    EXPECT_FALSE(edge1 < edge5);
    EXPECT_FALSE(edge6 < edge1);
    EXPECT_TRUE(edge1 < edge6);
    
    Edge edge7{1, 1, 10};
    Edge edge8{1, 3, 10};
    EXPECT_TRUE(edge7 < edge1);
    EXPECT_FALSE(edge1 < edge7);
    EXPECT_FALSE(edge8 < edge1);
    EXPECT_TRUE(edge1 < edge8);
}

TEST_F(EdgeTest, OperadorIgualdad) {
    EXPECT_TRUE(edge1 == edge3);
    EXPECT_FALSE(edge1 == edge2);
    EXPECT_FALSE(edge1 == edge4);

    Edge edge5{2, 2, 10};
    EXPECT_FALSE(edge1 == edge5);

    Edge edge6{1, 3, 10};
    EXPECT_FALSE(edge1 == edge6);

    Edge edge7{1, 2, 15};
    EXPECT_FALSE(edge1 == edge7);

    EXPECT_TRUE(edge3 == edge1);
}

TEST_F(EdgeTest, CreacionDeArista) {
    Edge e{5, 6, 15};
    EXPECT_EQ(e.node1, 5);
    EXPECT_EQ(e.node2, 6);
    EXPECT_EQ(e.weight, 15);

    Edge e1{0, 0, 0};
    EXPECT_EQ(e1.node1, 0);
    EXPECT_EQ(e1.node2, 0);
    EXPECT_EQ(e1.weight, 0);

    Edge e2{INT_MAX, INT_MAX, INT_MAX};
    EXPECT_EQ(e2.node1, INT_MAX);
    EXPECT_EQ(e2.node2, INT_MAX);
    EXPECT_EQ(e2.weight, INT_MAX);

    Edge e3{-1, -1, -1};
    EXPECT_EQ(e3.node1, -1);
    EXPECT_EQ(e3.node2, -1);
    EXPECT_EQ(e3.weight, -1);
}
