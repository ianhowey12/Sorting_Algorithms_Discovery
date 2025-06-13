#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <random>

using namespace std;

class Sort {
public:

    int index;
    string name;

    bool isComparison;
    string methodType;
    bool preservesOrder;
    bool inPlace;
    bool canSortDecimals;

    string bestTimeComplexity;
    string averageTimeComplexity;
    string worstTimeComplexity;

    string spaceComplexity;

    Sort() {}

    Sort(int index, string name, bool isComparison, string methodType, bool preservesOrder, bool inPlace, bool canSortDecimals, string bestTimeComplexity, string averageTimeComplexity, string worstTimeComplexity, string spaceComplexity);

    void printInfo();
};

struct Node;


class Sorter {
public:

    Sorter() {}

    void sortVector(int index, vector<int>& v, int n);

    void bubble(vector<int>& v, int n);

    void exchange(vector<int>& v, int n);

    void selection(vector<int>& v, int n);

    void double_selection(vector<int>& v, int n);

    void insertion(vector<int>& v, int n);

    int binary_search(vector<int>& v, int c, int l, int r);

    void binary_insertion(vector<int>& v, int n);

    void mer(vector<int>& v, int l, int m, int r);

    void merg(vector<int>& v, int l, int r);

    void merge(vector<int>& v, int n);

    void mer_3(vector<int>& v, int l, int m1, int m2, int r);

    void merg_3(vector<int>& v, int l, int r);

    void merge_3(vector<int>& v, int n);

    void mer_4(vector<int>& v, int l, int m1, int m2, int m3, int r);

    void merg_4(vector<int>& v, int l, int r);

    void merge_4(vector<int>& v, int n);

    int partition_lomuto(vector<int>& v, int l, int r);

    void qui_lomuto(vector<int>& v, int l, int r);

    void quick_lomuto(vector<int>& v, int n);

    int partition_hoare(vector<int>& v, int L, int R);

    void qui_hoare(vector<int>& v, int l, int r);

    void quick_hoare(vector<int>& v, int n);

    int partition_naive(vector<int>& v, int l, int r);

    void qui_naive(vector<int>& v, int l, int r);

    void quick_naive(vector<int>& v, int n);

    void heap_build(vector<int>& v, int n, int i);

    void heap(vector<int>& v, int n);

    void comb(vector<int>& v, int n);

    void shell(vector<int>& v, int n);

    void bogo(vector<int>& v, int n);

    void pigeonhole(vector<int>& v, int n);

    void counting(vector<int>& v, int n);

    void bucket_insertion(vector<int>& b);

    void bucket(vector<int>& v, int n, int N);

    void lsd_count(vector<int>& v, int n, int power, int b);

    void lsd(vector<int>& v, int n, int b);

    int numDigits(int n, int b);

    int pow(int b, int e);

    void msd_count(vector<int>& v, int left, int right, int digitPos, int b);

    int msd_power(vector<int>& v, int n, int b);

    void msd(vector<int>& v, int n, int b);

    void in_place_lsd(vector<int>& v, int n, int b);

    void in_place_msd_count(vector<int>& v, int l, int r, int digitPos, int b);

    int in_place_msd_power(vector<int>& v, int n, int b);

    void in_place_msd(vector<int>& v, int n, int b);

    void intro_insertion(vector<int>& v, int begin, int end);

    void intro_heap_build(vector<int>& v, int n, int i, int begin);

    void intro_heap(vector<int>& v, int begin, int end);

    int intro_partition(vector<int>& v, int begin, int end);

    void intro_sub(vector<int>& v, int l, int r, int limit);

    void intro(vector<int>& v, int n);

    void cycle(vector<int>& v, int n);

    void tim_insertion(vector<int>& v, int l, int r);

    void tim_merge(vector<int>& v, int l, int m, int r);

    void tim(vector<int>& v, int n, int b);

    void iterative_mer(vector<int>& v, int l, int m, int r);

    void iterative_merge(vector<int>& v, int n);

    void in_place_mer(vector<int>& v, int l, int m, int r);

    void in_place_merg(vector<int>& v, int l, int r);

    void in_place_merge(vector<int>& v, int n);

    void weave_swap(vector<int>& v, int i, int j, bool up);

    void weave_merge(vector<int>& v, int l, int s, bool up);

    void weav(vector<int>& v, int l, int s, bool up);

    void weave(vector<int>& v, int n);

    void iterative_weav(vector<int>& v, int l, int s, bool up);

    void iterative_weave(vector<int>& v, int n);

    void rotate_mer(vector<int>& v, int l, int m, int r);

    void rotate_merg(vector<int>& v, int l, int r);

    void rotate_merge(vector<int>& v, int n);

    void block_rotate(vector<int>& v, int l, int m, int r);

    int block_search(const vector<int>& v, int l, int r, int key);

    void block_merge(vector<int>& v, int l, int m, int r);

    void bloc(vector<int>& v, int l, int r);

    void block(vector<int>& v, int n);

    void iterative_block_merge(vector<int>& v, int l, int m, int r);

    void iterative_block(vector<int>& v, int n);

    void wiki(vector<int>& v, int n);

    void grail_sub(vector<int>& v);

    void grail(vector<int>& v, int n);

    void stooge_sub(vector<int>& v, int l, int r);

    void stooge(vector<int>& v, int n);

    void weak_heap_sift(vector<int>& v, int i, int n);
    
    void weak_heap_build(vector<int>& v, int n);

    void weak_heap(vector<int>& v, int n);
    
    int leonardo(int k);

    void build_heap(vector<int>& v, int l, int r);

    void smooth(vector<int>& v, int n);

    void poplar_heap_build(vector<int>& v, int n, int i);

    void poplar_heap(vector<int>& v, int n);

    int binary_partition(vector<int>& v, int l, int r);

    void binary_qui(vector<int>& v, int l, int r);

    void binary_quick(vector<int>& v, int n);

    void pancake_flip(vector<int>& v, int r);

    void pancake(vector<int>& v, int n);

    void cocktail(vector<int>& v, int n);

    void odd_even(vector<int>& v, int n);

    bool circle_pass(vector<int>& v, int l, int r);

    void circle(vector<int>& v, int n);

    void merge_insert_insertion(vector<int>& v, int p, int q);

    void merge_insert_merge(vector<int>& v, int l, int m, int r);

    void merge_insert(vector<int>& v, int l, int r, int b);

    void merge_insertion(vector<int>& v, int n, int b);

    Node* tree_insert(Node* root, int v);

    void tree_traverse(Node* root, vector<int>& sorted);

    void tree(vector<int>& v, int n);

    void tournament_build(vector<int>& tree, vector<int>& v, int n);

    void tournament_update(vector<int>& tree, int pos, int value, int n);

    int tournament_min(vector<int>& tree, int n);

    void tournament(vector<int>& v, int n);

    typedef pair<int, int> Pair;

    void stable_tournament_build(vector<Pair>& tree, const vector<Pair>& v, int n);

    void stable_tournament_update(vector<Pair>& tree, int pos, Pair value, int n);

    int stable_tournament_min(const vector<Pair>& tree, int n, Pair minVal);

    void stable_tournament(vector<int>& v, int n);

    void gnome(vector<int>& v, int n);

    void library(vector<int>& v, int n);

    vector<int> strand_merge(vector<int>& a, vector<int>& b);

    void strand(vector<int>& v, int n);

    void patience_merge(vector<int>& v, vector<vector<int>>& p);

    void patience(vector<int>& v, int n);

    void bitonic_merge(vector<int>& v, int low, int count, bool dir);
    
    void bitonic_sort(vector<int>& v, int low, int count, bool dir);

    void bitonic(vector<int>& v, int n);

    void odd_even_mer(vector<int>& v, int l, int n, int r);

    void odd_even_sort(vector<int>& v, int l, int n);

    void odd_even_merge(vector<int>& v, int n);

    void pairwise_network_merge(vector<int>& v, int l, int r);

    void pairwise_network_sort(vector<int>& v, int l, int r);

    void pairwise_network(vector<int>& v, int n);

    int spread_digits(int value);

    void spread_sort(vector<int>& v, int l, int r, int shift, int b);

    void spread(vector<int>& v, int n, int b);

    void flash(vector<int>& v, int n, float b);

    void bead(vector<int>& v, int n);

    void spaghetti(vector<int>& v, int n);

    void ska_digit(vector<int>& v, vector<int>& temp, int byte);

    void ska(vector<int>& v, int n);

};
