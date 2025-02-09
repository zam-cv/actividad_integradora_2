#include <gtest/gtest.h>
#include "DisjointSet.h"

TEST(disjoint_set_test, basic_operations) {
    DisjointSet ds(5);

    for (int i = 0; i < 5; i++) {
        EXPECT_EQ(ds.FindSet(i), i);
    }

    ds.UnionSets(0, 1);
    EXPECT_EQ(ds.FindSet(0), ds.FindSet(1));

    ds.UnionSets(2, 3);
    EXPECT_EQ(ds.FindSet(2), ds.FindSet(3));

    EXPECT_EQ(ds.FindSet(4), 4);

    ds.UnionSets(0, 2);
    EXPECT_EQ(ds.FindSet(0), ds.FindSet(2));
    EXPECT_EQ(ds.FindSet(1), ds.FindSet(3));

    int rep = ds.FindSet(0);
    EXPECT_EQ(ds.FindSet(1), rep);
    EXPECT_EQ(ds.FindSet(2), rep);
    EXPECT_EQ(ds.FindSet(3), rep);
}

TEST(disjoint_set_test, empty_set) {
    DisjointSet ds(0);
    EXPECT_EQ(ds.parent.size(), 0);
    EXPECT_EQ(ds.rank.size(), 0);
}

TEST(disjoint_set_test, union_by_rank) {
    DisjointSet ds(6);
    
    ds.UnionSets(0, 1);
    ds.UnionSets(1, 2);  
    
    ds.UnionSets(3, 4);  
    
    ds.UnionSets(2, 3);
    
    int rep = ds.FindSet(0);
    EXPECT_EQ(ds.FindSet(1), rep);
    EXPECT_EQ(ds.FindSet(2), rep);
    EXPECT_EQ(ds.FindSet(3), rep);
    EXPECT_EQ(ds.FindSet(4), rep);
    
    EXPECT_EQ(ds.FindSet(5), 5);
}

TEST(disjoint_set_test, path_compression) {
    DisjointSet ds(5);
    
    ds.parent[1] = 0;
    ds.parent[2] = 1;
    ds.parent[3] = 2;
    ds.parent[4] = 3;
    
    EXPECT_EQ(ds.FindSet(4), 0);
    
    EXPECT_EQ(ds.parent[4], 0);
    EXPECT_EQ(ds.parent[3], 0);
    EXPECT_EQ(ds.parent[2], 0);
    EXPECT_EQ(ds.parent[1], 0);
}

TEST(disjoint_set_test, unite_same_rank) {
    DisjointSet ds(4);
    
    ds.UnionSets(0, 1);
    ds.UnionSets(2, 3);
    
    ds.UnionSets(1, 2);
    
    int rep = ds.FindSet(0);
    EXPECT_EQ(ds.FindSet(1), rep);
    EXPECT_EQ(ds.FindSet(2), rep);
    EXPECT_EQ(ds.FindSet(3), rep);
    
    EXPECT_EQ(ds.rank[rep], 2);
}

TEST(disjoint_set_test, multiple_unions) {
    DisjointSet ds(8);
    
    ds.UnionSets(0, 1);
    ds.UnionSets(2, 3);
    ds.UnionSets(4, 5);
    ds.UnionSets(6, 7);
    
    ds.UnionSets(0, 2);
    ds.UnionSets(4, 6);
    
    ds.UnionSets(0, 4);
    
    int rep = ds.FindSet(0);
    for(int i = 0; i < 8; i++) {
        EXPECT_EQ(ds.FindSet(i), rep);
    }
}

TEST(disjoint_set_test, self_union) {
    DisjointSet ds(3);
    
    ds.UnionSets(1, 1);
    EXPECT_EQ(ds.FindSet(1), 1);
    EXPECT_EQ(ds.rank[1], 0);
}

TEST(disjoint_set_test, complex_path_compression) {
    DisjointSet ds(7);

    ds.parent[1] = 0;
    ds.parent[2] = 1;
    ds.parent[3] = 2;
    ds.parent[4] = 3;
    ds.parent[5] = 4;
    ds.parent[6] = 5;

    EXPECT_EQ(ds.FindSet(6), 0);

    for(int i = 1; i <= 6; i++) {
        EXPECT_EQ(ds.parent[i], 0);
    }
}
