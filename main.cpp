#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <chrono>

using namespace std;

#include "sort.hpp"

#define numSorts 82

vector<Sort> sorts(numSorts);

// vector to be sorted and its size
vector<int> v(0, 0);
int n = 0;

Sort::Sort(int index, string name, bool isComparison, string methodType, bool preservesOrder, bool inPlace, bool canSortDecimals, string bestTimeComplexity, string averageTimeComplexity, string worstTimeComplexity, string spaceComplexity) {
    this->index = index;
    this->name = name;
    this->isComparison = isComparison;
    this->methodType = methodType;
    this->preservesOrder = preservesOrder;
    this->inPlace = inPlace;
    this->canSortDecimals = canSortDecimals;
    this->bestTimeComplexity = bestTimeComplexity;
    this->averageTimeComplexity = averageTimeComplexity;
    this->worstTimeComplexity = worstTimeComplexity;
    this->spaceComplexity = spaceComplexity;
}

void Sort::printInfo() {
    cout << "Sorting Algorithm #" << index << ": " << name << "\nIs a comparison sort: " << isComparison << "\nAlgorithm Type: " << methodType << "\nPreserves order of duplicate elements (stable): "
    << preservesOrder << "\nSorts elements in-place (without additional storage): " << inPlace << "\nCan be modified to sort decimals: " << canSortDecimals << "\nBest-case time complexity: " << bestTimeComplexity
    << "\nAverage-case time complexity: " << averageTimeComplexity << "\nWorst-case time complexity: " << worstTimeComplexity << "\nSpace complexity: " << spaceComplexity << "\n\n";
}

void Sorter::sortVector(int index, vector<int>& v, int n) {
    switch (index) {
    case 0: // bubble
        bubble(v, n);
        break;
    case 1: // exchange
        exchange(v, n);
        break;
    case 2: // selection
        selection(v, n);
        break;
    case 3: // double selection
        double_selection(v, n);
        break;
    case 4: // insertion
        insertion(v, n);
        break;
    case 5: // binary insertion
        binary_insertion(v, n);
        break;
    case 6: // merge
        merge(v, n);
        break;
    case 7: // three-way merge
        merge_3(v, n);
        break;
    case 8: // four-way merge
        merge_4(v, n);
        break;
    case 9: // quick, lomuto partition
        quick_lomuto(v, n);
        break;
    case 10: // quick, hoare partition
        quick_hoare(v, n);
        break;
    case 11: // quick, naive partition
        quick_naive(v, n);
        break;
    case 12: // heap
        heap(v, n);
        break;
    case 13: // comb
        comb(v, n);
        break;
    case 14: // shell
        shell(v, n);
        break;
    case 15: // bogo
        bogo(v, n);
        for(int j=0;j<n;j++){
        cout << v[j] << " ";
        }
        cout << "\n";
        break;

    case 16: // pigeonhole
        pigeonhole(v, n);
        break;
    case 17: // counting
        counting(v, n);
        break;
    case 18: // bucket, 4 buckets
        bucket(v, n, 4);
        break;
    case 19: // bucket, 16 buckets
        bucket(v, n, 16);
        break;
    case 20: // bucket, 64 buckets
        bucket(v, n, 64);
        break;
    case 21: // bucket, N buckets
        bucket(v, n, n);
        break;
    case 22: // lsd radix, base 2
        lsd(v, n, 2);
        break;
    case 23: // lsd radix, base 4
        lsd(v, n, 4);
        break;
    case 24: // lsd radix, base 8
        lsd(v, n, 8);
        break;
    case 25: // lsd radix, base 16
        lsd(v, n, 16);
        break;
    case 26: // msd radix, base 2
        msd(v, n, 2);
        break;
    case 27: // msd radix, base 4
        msd(v, n, 4);
        break;
    case 28: // msd radix, base 8
        msd(v, n, 8);
        break;
    case 29: // msd radix, base 16
        msd(v, n, 16);
        break;
    case 30: // in place lsd radix, base 2
        in_place_lsd(v, n, 2);
        break;
    case 31: // in place lsd radix, base 4
        in_place_lsd(v, n, 4);
        break;
    case 32: // in place lsd radix, base 8
        in_place_lsd(v, n, 8);
        break;
    case 33: // in place lsd radix, base 16
        in_place_lsd(v, n, 16);
        break;
    case 34: // in place msd radix, base 2
        in_place_msd(v, n, 2);
        break;
    case 35: // in place msd radix, base 4
        in_place_msd(v, n, 4);
        break;
    case 36: // in place msd radix, base 8
        in_place_msd(v, n, 8);
        break;
    case 37: // in place msd radix, base 16
        in_place_msd(v, n, 16);
        break;

    case 38: // intro
        intro(v, n);
        break;
    case 39: // cycle
        cycle(v, n);
        break;
    case 40: // tim, group size 4
        tim(v, n, 4);
        break;
    case 41: // tim, group size 16
        tim(v, n, 16);
        break;
    case 42: // tim, group size 64
        tim(v, n, 64);
        break;
    case 43: // iterative merge
        iterative_merge(v, n);
        break;
    case 44: // naive in-place merge
        in_place_merge(v, n);
        break;
    case 45: // weave
        weave(v, n);
        break;
    case 46: // iterative weave
        iterative_weave(v, n);
        break;
    case 47: // rotate merge
        rotate_merge(v, n);
        break;
    case 48: // block
        block(v, n);
        break;
    case 49: // iterative block
        iterative_block(v, n);
        break;
    case 50: // wiki
        wiki(v, n);
        break;
    case 51: // grail
        grail(v, n);
        break;
    case 52: // stooge
        stooge(v, n);
        break;
    case 53: // weak heap
        weak_heap(v, n);
        break;
    case 54: // smooth
        smooth(v, n);
        break;
    case 55: // poplar heap
        poplar_heap(v, n);
        break;
    case 56: // binary quick
        binary_quick(v, n);
        break;
    case 57: // pancake
        pancake(v, n);
        break;
    case 58: // cocktail shaker
        cocktail(v, n);
        break;
    case 59: // odd-even
        odd_even(v, n);
        break;
    case 60: // circle
        circle(v, n);
        break;
    case 61: // merge-insertion size 8
        merge_insertion(v, n, 8);
        break;
    case 62: // tree
        tree(v, n);
        break;
    case 63: // tournament
        tournament(v, n);
        break;
    case 64: // stable tournament
        stable_tournament(v, n);
        break;
    case 65: // gnome
        gnome(v, n);
        break;
    case 66: // library
        library(v, n);
        break;
    case 67: // strand
        strand(v, n);
        break;
    case 68: // patience
        patience(v, n);
        break;
    case 69: // bitonic
        bitonic(v, n);
        break;
    case 70: // odd-even merge
        odd_even_merge(v, n);
        break;
    case 71: // pairwise network
        pairwise_network(v, n);
        break;

    case 72: // spread, 1 bit base
        spread(v, n, 1);
        break;
    case 73: // spread, 2 bit base
        spread(v, n, 2);
        break;
    case 74: // spread, 4 bit base
        spread(v, n, 4);
        break;
    case 75: // spread, 8 bit base
        spread(v, n, 8);
        break;
    case 76: // flash, 0.3n buckets
        flash(v, n, 0.3);
        break;
    case 77: // flash, 0.4n buckets
        flash(v, n, 0.4);
        break;
    case 78: // flash, 0.5n buckets
        flash(v, n, 0.5);
        break;
    case 79: // bead (gravity)
        bead(v, n);
        break;
    case 80: // spaghetti
        spaghetti(v, n);
        break;
    case 81: // ska
        ska(v, n);
        break;
    }
}

void Sorter::bubble(vector<int>& v, int n) {
    bool swapped;

    for (int i = 0; i < n - 1; i++) {
        swapped = 0;
        int end = n - i - 1;
        for (int j = 0; j < end; j++) {
            if (v[j] > v[j + 1]) {
                swap(v[j], v[j + 1]);
                swapped = 1;
            }
        }

        if (!swapped) break;
    }
}

void Sorter::exchange(vector<int>& v, int n){
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (v[j] < v[i]) {
                int temp = v[i];
                v[i] = v[j];
                v[j] = temp;
            }
        }
    }
}

void Sorter::selection(vector<int>& v, int n) {

    for (int i = 0; i < n; i++) {
        int minI = i;

        for (int j = i + 1; j < n; j++) {
            if (v[j] < v[minI]) {
                minI = j;
            }
        }

        swap(v[i], v[minI]);
    }
}

void Sorter::double_selection(vector<int>& v, int n)
{
    for (int i = 0, j = n - 1; i < j; i++, j--) {
        int mini = v[i], maxi = v[i], minI = i, maxI = i;

        for (int k = i; k <= j; k++) {
            if (v[k] > maxi) {
                maxi = v[k];
                maxI = k;
            }
            else if (v[k] < mini) {
                mini = v[k];
                minI = k;
            }
        }

        // shifting the min
        swap(v[i], v[minI]);

        // shifting the max. Same condition happens if we shifted the max to arr[minI] in the previous swap
        if (v[minI] == maxi)
            swap(v[j], v[minI]);
        else
            swap(v[j], v[maxI]);
    }
}

void Sorter::insertion(vector<int>& v, int n){
    for (int i = 1; i < n; i++) {
        int val = v[i];
        int j = i - 1;

        while (j >= 0 && v[j] > val) {
            v[j + 1] = v[j];
            j--;
        }
        v[j + 1] = val;
    }
}

int Sorter::binary_search(vector<int>& v, int c, int l, int r){
    if (r <= l) return (c > v[l]) ? (l + 1) : l;

    int m = (l + r) / 2;

    if (c == v[m]) return m + 1;

    if (c > v[m]) return binary_search(v, c, m + 1, r);
    return binary_search(v, c, l, m - 1);
}

void Sorter::binary_insertion(vector<int>& v, int n){
    int i, loc, j, c;

    for (i = 1; i < n; ++i){
        j = i - 1;
        c = v[i];

        // find location where c should be inserted
        loc = binary_search(v, c, 0, j);

        // Move all elements after location to create space
        while (j >= loc)
        {
            v[j + 1] = v[j];
            j--;
        }
        v[j + 1] = c;
    }
}

void Sorter::mer(vector<int>&v, int l, int m, int r){
    int n1 = m - l + 1;
    int n2 = r - m;

    // Create temp vectors
    vector<int> L(n1), R(n2);

    // Copy data to temp vectors L and R
    for (int i = 0; i < n1; i++) L[i] = v[l + i];
    for (int j = 0; j < n2; j++) R[j] = v[m + 1 + j];

    int i = 0, j = 0;
    int k = l;

    // Merge the temp vectors back into v[left..right]
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            v[k] = L[i];
            i++;
        }
        else {
            v[k] = R[j];
            j++;
        }
        k++;
    }

    // Copy the remaining elements of L
    while (i < n1) {
        v[k] = L[i];
        i++;
        k++;
    }

    // Copy the remaining elements of R
    while (j < n2) {
        v[k] = R[j];
        j++;
        k++;
    }
}

void Sorter::merg(vector<int>& v, int l, int r) {
    if (l >= r) return;

    int m = l + (r - l) / 2;
    merg(v, l, m);
    merg(v, m + 1, r);
    mer(v, l, m, r);
}

void Sorter::merge(vector<int>& v, int n){
    merg(v, 0, n);
}

void Sorter::mer_3(vector<int>& v, int l, int m1, int m2, int r) {

    // Subarrays
    int n1 = m1 - l + 1;
    int n2 = m2 - m1;
    int n3 = r - m2;
    vector<int> v1(n1), v2(n2), v3(n3);

    // Copy data to temporary arrays
    for (int i = 0; i < n1; i++) {
        v1[i] = v[l + i];
    }
    for (int i = 0; i < n2; i++) {
        v2[i] = v[m1 + 1 + i];
    }
    for (int i = 0; i < n3; i++) {
        v3[i] = v[m2 + 1 + i];
    }

    // Merge sorted subarrays
    int i = 0, j = 0, k = 0, spot = l;
    while (i < n1 || j < n2 || k < n3) {
        int minVal = INT_MAX;
        int m = 0;

        // Find smallest
        if(i < n1){
            m = 0; minVal = v1[i];
        }
        if(j < n2 && v2[j] < minVal){
            m = 1; minVal = v2[j];
        }
        if(k < n3 && v3[k] < minVal){
            m = 2;
        }
        // Place and increment
        switch(m){
            case 0:
                v[spot++] = v1[i++];
                break;
            case 1:
                v[spot++] = v2[j++];
                break;
            case 2:
                v[spot++] = v3[k++];
                break;
        }
    }
}

void Sorter::merg_3(vector<int>& v, int l, int r) {

    if (l >= r) {
        return;
    }

    int m1 = l + (r - l) / 3;
    int m2 = l + 2 * (r - l) / 3;

    merg_3(v, l, m1);
    merg_3(v, m1 + 1, m2);
    merg_3(v, m2 + 1, r);

    mer_3(v, l, m1, m2, r);
}

void Sorter::merge_3(vector<int>& v, int n) {
    merg_3(v, 0, n - 1);
}

void Sorter::mer_4(vector<int>& v, int l, int m1, int m2, int m3, int r) {
    
    // Subarrays
    int n1 = m1 - l;
    int n2 = m2 - m1;
    int n3 = m3 - m2;
    int n4 = r - m3;
    vector<int> v1(n1), v2(n2), v3(n3), v4(n4);

    // Copy the subarrays into temporary vectors
    for (int i = 0; i < n1; i++) v1[i] = v[l + i];
    for (int i = 0; i < n2; i++) v2[i] = v[m1 + i];
    for (int i = 0; i < n3; i++) v3[i] = v[m2 + i];
    for (int i = 0; i < n4; i++) v4[i] = v[m3 + i];

    // Merge sorted subarrays
    int h = 0, i = 0, j = 0, k = 0, spot = l;
    while (h < n1 || i < n2 || j < n3 || k < n4) {
        int minVal = INT_MAX;
        int m = 0;

        // Find smallest
        if(h < n1){
            m = 0; minVal = v1[h];
        }
        if(i < n2 && v2[i] < minVal){
            m = 1; minVal = v2[i];
        }
        if(j < n3 && v3[j] < minVal){
            m = 2; minVal = v3[j];
        }
        if(k < n4 && v4[k] < minVal){
            m = 3;
        }
        // Place and increment
        switch(m){
            case 0:
                v[spot++] = v1[h++];
                break;
            case 1:
                v[spot++] = v2[i++];
                break;
            case 2:
                v[spot++] = v3[j++];
                break;
            case 3:
                v[spot++] = v4[k++];
                break;
        }
    }
}

void Sorter::merg_4(vector<int>& v, int l, int r) {
    if (r - l <= 1) return;

    int s = r - l;
    int m1 = l + s / 4;
    int m2 = l + s / 2;
    int m3 = l + 3 * s / 4;

    merg_4(v, l, m1);
    merg_4(v, m1, m2);
    merg_4(v, m2, m3);
    merg_4(v, m3, r);

    mer_4(v, l, m1, m2, m3, r);
}

void Sorter::merge_4(vector<int>& v, int n) {
    merg_4(v, 0, n);
}

int Sorter::partition_lomuto(vector<int>& v, int l, int r) {

    int pivot = v[r];

    // Index of smaller element and indicates the right position of pivot found so far
    int i = l - 1;

    // Traverse v[low..high] and move all smaller elements on left side. Elements from l to i are smaller after every iteration
    for (int j = l; j <= r - 1; j++) {
        if (v[j] < pivot) {
            i++;
            swap(v[i], v[j]);
        }
    }

    // Move pivot after smaller elements and return its position
    swap(v[i + 1], v[r]);
    return i + 1;
}

void Sorter::qui_lomuto(vector<int>& v, int l, int r) {
    if (l < r) {
        int i = partition_lomuto(v, l, r);
        qui_lomuto(v, l, i - 1);
        qui_lomuto(v, i + 1, r);
    }
}

void Sorter::quick_lomuto(vector<int>& v, int n) {
    qui_lomuto(v, 0, n-1);
}

int Sorter::partition_hoare(vector<int>& v, int L, int R) {
    int pivot = v[L];

    int l = L - 1;
    int r = R + 1;
 
    while (1) {
        l++;
        while (v[l] < pivot) {
            l++;
        }
 
        r--;
        while (v[r] > pivot) {
            r--;
        }

        if (l >= r) return r;

        swap(v[l], v[r]);
    }
}

void Sorter::qui_hoare(vector<int>& v, int l, int r) {
    if (l < r) { 
        int i = partition_hoare(v, l, r);
        qui_hoare(v, l, i); 
        qui_hoare(v, i + 1, r);
    }
}

void Sorter::quick_hoare(vector<int>& v, int n) {
    qui_hoare(v, 0, n-1);
}

int Sorter::partition_naive(vector<int>& v, int l, int r) {
    int n = r - l + 1;

    // Last element will be the pivot value
    int pivot = v[r];
    int p = 0;

    // create a temp array to store the elements in order
    vector<int> temp(n);
    int spot = 0;

    // Fill elements <= pivot
    for (int i = l; i < r; i++) {
        if (v[i] <= pivot) temp[spot++] = v[i];
    }
    // Move the pivot and record its new index within temp for returning
    p = l + spot;
    temp[spot++] = pivot;

    // Fill elements > pivot
    for (int i = l; i < r; i++) {
        if (v[i] > pivot) temp[spot++] = v[i];
    }

    // Copy elements from temp to the subset [i, j] of v
    for (int i = 0; i < n; i++) {
        v[l+i] = temp[i];
    }

    // Edge case
    if(p == r) return p - 1;
    return p;
}

void Sorter::qui_naive(vector<int>& v, int l, int r) {
    if (l < r) {
        int p = partition_naive(v, l, r);
        qui_naive(v, l, p);
        qui_naive(v, p + 1, r);
    }
}

void Sorter::quick_naive(vector<int>& v, int n) {
    qui_naive(v, 0, n-1);
}

void Sorter::heap_build(vector<int>& v, int n, int i) {

    // Initialize max as root
    int max = i;

    // children
    int l = 2 * i + 1;
    int r = 2 * i + 2;

    // Update max with children
    if (l < n && v[l] > v[max]) max = l;
    if (r < n && v[r] > v[max]) max = r;

    // If max is not root
    if (max != i) {
        swap(v[i], v[max]);

        // Make subset into sub-heap
        heap_build(v, n, max);
    }
}

void Sorter::heap(vector<int>& v, int n) {

    // Build heap (rearrange vector)
    for (int i = n / 2 - 1; i >= 0; i--)
        heap_build(v, n, i);

    // One by one extract an element from heap
    for (int i = n - 1; i > 0; i--) {

        // Move current root to end
        swap(v[0], v[i]);

        heap_build(v, i, 0);
    }
}

void Sorter::comb(vector<int>& v, int n)
{
    // Initialize gap
    int g = n;

    // Initialize swapped as true to make sure that loop runs
    bool swapped = 1;

    while (g > 1 || swapped){ // IDK IF SWAPPED SHOULD BE CHANGED TO !SWAPPED HERE
        // Find next gap
        g = (g * 10) / 13;
        // rule of 11
        if(g == 9 || g == 10) g = 11;
        if (g < 1) g = 1;

        // Initialize swapped as false so that we can check if swap happened or not
        swapped = 0;

        // Compare all elements with current gap
        for (int i = 0; i < n - g; i++){
            if (v[i] > v[i + g]){
                swap(v[i], v[i + g]);
                swapped = 1;
            }
        }
    }
}

void Sorter::shell(vector<int>& v, int n) {
    // Start with a large gap, then reduce it
    for (int g = n / 2; g > 0; g /= 2) {
        // Perform a gapped insertion sort for this gap size
        for (int i = g; i < n; i++) {
            int temp = v[i];
            int j;
            for (j = i; j >= g && v[j - g] > temp; j -= g) {
                v[j] = v[j - g];
            }
            v[j] = temp;
        }
    }
}

void Sorter::bogo(vector<int>& v, int n) {
    int N; bool sorted;
    while (1) {
        N = n;
        sorted = 1;

        // if array is not sorted then shuffle the array again
        while (--N > 0)
            if (v[N] < v[N - 1]) {
                sorted = 0;
                break;
            }
        if (sorted) return;

        for (int i = 0; i < n; i++) swap(v[i], v[rand() % n]);
    }
}

void Sorter::pigeonhole(vector<int>& v, int n) {
    int min = v[0], max = v[0];
    for (int i = 1; i < n; i++)
    {
        if (v[i] < min) min = v[i];
        if (v[i] > max) max = v[i];
    }
    int r = max - min + 1;

    vector<vector<int>> holes(r);

    // Put every element in its hole
    for (int i = 0; i < n; i++) holes[v[i] - min].push_back(v[i]);

    int index = 0;
    for (int i = 0; i < r; i++)
    {
        vector<int>::iterator it;
        for (it = holes[i].begin(); it != holes[i].end(); ++it) v[index++] = *it;
    }
}

void Sorter::counting(vector<int>& v, int n) {
    // WARNING: MAX ELEMENT CANNOT BE HUGE AND NO ELEMENTS CAN BE NEGATIVE because of c[v[i]] access.

    // Finding the maximum element of v
    int m = 0;

    for (int i = 0; i < n; i++)
        m = max(m, v[i]);

    // Initializing c with 0
    vector<int> c(m + 1, 0);

    // Mapping each element of v as an index of c
    for (int i = 0; i < n; i++)
        c[v[i]]++;

    // Calculating prefix sum at every index of c
    for (int i = 1; i <= m; i++)
        c[i] += c[i - 1];

    // Creating out from c
    vector<int> out(n);

    for (int i = n - 1; i >= 0; i--) {
        out[c[v[i]] - 1] = v[i];
        c[v[i]]--;
    }

    v = out;
}

void Sorter::bucket_insertion(vector<int>& b) {
    int n = size(b);
    for (int i = 1; i < n; i++) {
        int val = b[i];
        int j = i - 1;
        while (j >= 0 && b[j] > val) {
            b[j + 1] = b[j];
            j--;
        }
        b[j + 1] = val;
    }
}

void Sorter::bucket(vector<int>& v, int n, int N) {
    // Create N empty buckets
    vector<vector<int>> b(N);

    // Find min and max to determine bucket divisions
    int min = v[0];
    for (int i = 1; i < n; i++)
        if (v[i] < min) min = v[i];
    int max = v[0];
    for (int i = 1; i < n; i++)
        if (v[i] > max) max = v[i];
    int range = max - min;

    // Put array elements in buckets using their values
    for (int i = 0; i < n; i++) {
        int bi = 0.999 * (double)N * (double)(v[i] - min) / (double)range;
        b[bi].push_back(v[i]);
    }

    // Sort individual buckets using insertion sort
    for (int i = 0; i < N; i++) {
        bucket_insertion(b[i]);
    }

    // Concatenate all buckets into v
    int index = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < b[i].size(); j++) {
            v[index++] = b[i][j];
        }
    }
}

void Sorter::lsd_count(vector<int>& v, int n, int power, int b){
    vector<int> out(n);
    int i;
    vector<int> c(b, 0);

    // Store count of occurrences in c
    for (i = 0; i < n; i++)
        c[(v[i] / power) % b]++;

    // Get cumulative counts
    for (i = 1; i < b; i++)
        c[i] += c[i - 1];

    // Build the output array
    for (i = n - 1; i >= 0; i--) {
        int digit = (v[i] / power) % b;
        out[c[digit] - 1] = v[i];
        c[digit]--;
    }

    v = out;
}

void Sorter::lsd(vector<int>& v, int n, int b){
    // Find the maximum number to know number of digits
    int m = v[0];
    for (int i = 1; i < n; i++)
        if (v[i] > m) m = v[i];

    // Do counting sort for every digit's power
    for (int power = 1; m / power > 0; power *= b) lsd_count(v, n, power, b);
}

int Sorter::numDigits(int n, int b){
    int count = 0;
    while (n > 0){
        n /= b;
        count++;
    }
    return count;
}

int Sorter::pow(int b, int e) {
    if (e < 1) return 1;
    int n = b;
    for (int i = 1; i < e; i++) {
        n *= b;
    }
    return n;
}

void Sorter::msd_count(vector<int>& v, int l, int r, int power, int b) {
    if (l >= r || power == 0) return;

    vector<int> count(b+1, 0);

    // Count cumulatively
    for (int i = l; i <= r; i++) {
        int digit = (v[i] / power) % b;
        count[digit]++;
    }
    int temp = count[0];
    int temp2 = 0;
    count[0] = l;
    for (int i = 1; i <= b; i++) {
        temp2 = temp;
        temp = count[i];
        count[i] = temp2 + count[i - 1];
    }

    vector<int> o(r-l+1, 0);

    vector<int> index(b+1, 0);
    for (int i = 0; i <= b; i++) {
        index[i] = count[i];
    }

    // Rearrangement
    for (int i = l; i <= r; i++) {
        int digit = (v[i] / power) % b;
        o[count[digit] - l] = v[i];
        count[digit]++;
    }
    for (int i = l; i <= r; i++) {
        v[i] = o[i-l];
    }

    // Recurse on each bucket
    for (int i = 1; i <= b; i++) {
        msd_count(v, index[i-1], index[i] - 1, power / b, b);
    }
}

int Sorter::msd_power(vector<int>& v, int n, int b) {
    int maxVal = v[0];
    for(int i=1;i<n;i++){
        maxVal = max(maxVal, v[i]);
    }
    int power = 1;
    while (maxVal >= b) {
        power *= b;
        maxVal /= b;
    }
    return power;
}

void Sorter::msd(vector<int>& v, int n, int b) {
    if (n <= 1) return;
    msd_count(v, 0, n - 1, msd_power(v, n, b), b);
}

void Sorter::in_place_lsd(vector<int>& v, int n, int b) {
    if (n == 0) return;

    // Find max number to determine number of digits
    int maxVal = v[0];
    for(int i=1;i<n;i++){
        maxVal = max(maxVal, v[i]);
    }

    // Temporary buffer
    vector<int> output(n);

    int exp = 1;

    while (maxVal / exp > 0) {
        vector<int> count(b, 0);

        // Count occurrences of each digit
        for (int i = 0; i < n; ++i) {
            int digit = (v[i] / exp) % b;
            count[digit]++;
        }

        // Accumulate counts
        for (int i = 1; i < b; ++i) {
            count[i] += count[i - 1];
        }

        // Build output array (stable sort)
        for (int i = n - 1; i >= 0; --i) {
            int digit = (v[i] / exp) % b;
            output[--count[digit]] = v[i];
        }

        // Copy output back to original vector
        for (int i = 0; i < n; ++i) {
            v[i] = output[i];
        }

        exp *= b;
    }
}

void Sorter::in_place_msd_count(vector<int>& v, int l, int r, int power, int b) {
    if (l >= r || power == 0) return;

    vector<int> count(b, 0);
    vector<int> index(b, 0);

    // Count cumulatively
    for (int i = l; i <= r; i++) {
        int digit = (v[i] / power) % b;
        count[digit]++;
    }
    index[0] = l;
    for (int i = 1; i < b; i++) {
        index[i] = index[i - 1] + count[i - 1];
    }

    // In-place rearrangement using cycle leader iteration
    for (int i = 0; i < b; i++) {
        while (count[i] > 0) {
            int from = index[i];
            int val = v[from];
            int digit = (val / power) % b;

            while (digit != i) {
                swap(val, v[index[digit]]);
                index[digit]++;
                count[digit]--;
                digit = (val / power) % b;
            }

            v[from] = val;
            index[i]++;
            count[i]--;
        }
    }

    // Recurse on each bucket
    int start = l;
    for (int i = 0; i < b; ++i) {
        int size = (i == 0 ? index[0] : index[i]) - start;
        if (size > 1) {
            in_place_msd_count(v, start, index[i] - 1, power / b, b);
        }
        start = index[i];
    }
}

int Sorter::in_place_msd_power(vector<int>& v, int n, int b) {
    int maxVal = v[0];
    for(int i=1;i<n;i++){
        maxVal = max(maxVal, v[i]);
    }
    int power = 1;
    while (maxVal >= b) {
        power *= b;
        maxVal /= b;
    }
    return power;
}

void Sorter::in_place_msd(vector<int>& v, int n, int b) {
    if (n <= 1) return;
    in_place_msd_count(v, 0, n - 1, in_place_msd_power(v, n, b), b);
}

void Sorter::intro_insertion(vector<int>& v, int begin, int end) {
    for (int i = begin + 1; i < end; ++i) {
        int key = v[i];
        int j = i - 1;
        while (j >= begin && v[j] > key) {
            v[j + 1] = v[j];
            j--;
        }
        v[j + 1] = key;
    }
}

void Sorter::intro_heap_build(vector<int>& v, int n, int i, int begin) {
    int largest = i;
    int l = 2 * (i - begin) + 1 + begin;
    int r = 2 * (i - begin) + 2 + begin;

    if (l < n && v[l] > v[largest])
        largest = l;
    if (r < n && v[r] > v[largest])
        largest = r;
    if (largest != i) {
        swap(v[i], v[largest]);
        intro_heap_build(v, n, largest, begin);
    }
}

void Sorter::intro_heap(vector<int>& v, int begin, int end) {
    int n = end;
    // Build heap
    for (int i = (begin + (n - begin) / 2) - 1; i >= begin; --i)
        intro_heap_build(v, n, i, begin);
    // Extract elements
    for (int i = n - 1; i > begin; --i) {
        swap(v[begin], v[i]);
        intro_heap_build(v, i, begin, begin);
    }
}

int Sorter::intro_partition(vector<int>& v, int begin, int end) {
    int pivot = v[end - 1];
    int i = begin - 1;
    for (int j = begin; j < end - 1; ++j) {
        if (v[j] <= pivot) {
            ++i;
            swap(v[i], v[j]);
        }
    }
    swap(v[i + 1], v[end - 1]);
    return i + 1;
}

void Sorter::intro_sub(vector<int>& v, int l, int r, int limit) {
    int size = r - l;

    if (size < 16) {
        intro_insertion(v, l, r);
        return;
    }

    if (limit == 0) {
        intro_heap(v, l, r);
        return;
    }

    int pivot = intro_partition(v, l, r);
    intro_sub(v, l, pivot, limit - 1);
    intro_sub(v, pivot + 1, r, limit - 1);
}

void Sorter::intro(vector<int>& v, int n) {
    int limit = 2 * log(n);
    intro_sub(v, 0, n, limit);
}

void Sorter::cycle(vector<int>& v, int n) {
    // count number of memory writes
    int writes = 0;

    for (int i = 0; i <= n - 2; i++) {
        // initialize item as starting point
        int item = v[i];

        // find position of item (count all smaller elements on right side of item)
        int pos = i;
        for (int j = i + 1; j < n; j++)
            if (v[j] < item) pos++;

        // if item is already in correct position
        if (pos == i) continue;

        // ignore all duplicate elements
        while (item == v[pos]) pos += 1;

        // put the item in its right position
        if (pos != i) {
            swap(item, v[pos]);
            writes++;
        }

        // rotate rest of the cycle
        while (pos != i) {
            pos = i;

            // find position where we put the element
            for (int j = i + 1; j < n; j++)
                if (v[j] < item)
                    pos += 1;

            // ignore all duplicate elements
            while (item == v[pos])
                pos += 1;

            // put the item in its right position
            if (item != v[pos]) {
                swap(item, v[pos]);
                writes++;
            }
        }
    }

    // Number of memory writes or swaps
    // cout << writes << endl;
}

void Sorter::tim_insertion(vector<int>& v, int l, int r){
    for (int i = l + 1; i <= r; i++) {
        int temp = v[i];
        int j = i - 1;
        while (j >= l && v[j] > temp) {
            v[j + 1] = v[j];
            j--;
        }
        v[j + 1] = temp;
    }
}

void Sorter::tim_merge(vector<int>& v, int l, int m, int r){
    // Original array is broken in two
    int s1 = m - l + 1, s2 = r - m;
    vector<int> L(s1), R(s2);
    for (int i = 0; i < s1; i++)
        L[i] = v[l + i];
    for (int i = 0; i < s2; i++)
        R[i] = v[m + 1 + i];

    int i = 0;
    int j = 0;
    int k = l;

    // After comparing, we merge those two array in larger sub array 
    while (i < s1 && j < s2) {
        if (L[i] <= R[j]) {
            v[k] = L[i];
            i++;
        }
        else {
            v[k] = R[j];
            j++;
        }
        k++;
    }

    // Copy remaining elements of left
    while (i < s1) {
        v[k] = L[i];
        k++;
        i++;
    }

    // Copy remaining elements of right
    while (j < s2) {
        v[k] = R[j];
        k++;
        j++;
    }
}

void Sorter::tim(vector<int>& v, int n, int b){
    // Sort individual subarrays of size b or less
    for (int i = 0; i < n; i += b)
        tim_insertion(v, i, min((i + b - 1), (n - 1)));

    // Start merging from size b
    for (int s = b; s < n; s *= 2) {

        for (int l = 0; l < n; l += 2 * s) {

            int m = l + s - 1;
            int r = min((l + 2 * s - 1), (n - 1));

            if (m < r) tim_merge(v, l, m, r);
        }
    }
}

void Sorter::iterative_mer(vector<int>& v, int l, int m, int r) {
    vector<int> temp(r - l + 1);

    int i = l;
    int j = m + 1;
    int k = 0;

    // Merge both halves into temp
    while (i <= m && j <= r) {
        if (v[i] <= v[j])
            temp[k++] = v[i++];
        else
            temp[k++] = v[j++];
    }

    // Copy remaining elements from left half
    while (i <= m)
        temp[k++] = v[i++];

    // Copy remaining elements from right half
    while (j <= r)
        temp[k++] = v[j++];

    // Copy sorted elements back to original vector
    for (int i = 0; i < temp.size(); ++i)
        v[l + i] = temp[i];
}

void Sorter::iterative_merge(vector<int>& v, int n) {

    // w is size of subarrays to merge, starts at 1 and doubles each iteration
    for (int w = 1; w < n; w *= 2) {
        for (int l = 0; l < n; l += 2 * w) {
            int m = min(l + w - 1, n - 1);
            int r = min(l + 2 * w - 1, n - 1);

            if (m < r) {
                iterative_mer(v, l, m, r);
            }
        }
    }
}

void Sorter::in_place_mer(vector<int>& v, int l, int m, int r) {
    int i = l;
    int j = m + 1;

    while (i <= m && j <= r) {
        if (v[i] <= v[j]) {
            i++;
        }
        else {
            int value = v[j];
            int index = j;

            // Shift all elements between i and j one step to the right
            while (index > i) {
                v[index] = v[index - 1];
                index--;
            }

            v[i] = value;

            // Update all pointers
            i++;
            m++;
            j++;
        }
    }
}

void Sorter::in_place_merg(vector<int>& v, int l, int r) {
    if (l >= r) return;

    int m = l + (r - l) / 2;
    in_place_merg(v, l, m);
    in_place_merg(v, m + 1, r);
    in_place_mer(v, l, m, r);
}

void Sorter::in_place_merge(vector<int>& v, int n) {
    in_place_merg(v, 0, n - 1);
}

void Sorter::weave_swap(vector<int>& v, int i, int j, bool up) {
    if (up) {
        if (v[i] > v[j]) swap(v[i], v[j]);
    }
    else {
        if (v[i] < v[j]) swap(v[i], v[j]);
    }
}

void Sorter::weave_merge(vector<int>& v, int l, int s, bool up) {
    if (s == 1) return;

    int m = s / 2;
    for (int i = 0; i < m; ++i) {
        weave_swap(v, l + i, l + m + i, up);
    }

    weave_merge(v, l, m, up);
    weave_merge(v, l + m, m, up);
}

void Sorter::weav(vector<int>& v, int l, int s, bool up) {
    if (s <= 1) return;

    int m = s / 2;
    weav(v, l, m, 1);            // sort first half ascending
    weav(v, l + m, m, 0);     // sort second half descending

    weave_merge(v, l, s, up);     // merge whole section
}

void Sorter::weave(vector<int>& v, int n) {
    int power = 1;
    while (power < n) power *= 2;

    while (size(v) < power) {
        v.push_back(INT_MAX);
    }

    weav(v, 0, size(v), 1);

    // Remove padded values
    while (size(v) > n) {
        v.pop_back();
    }
}

void Sorter::iterative_weav(vector<int>& v, int l, int s, bool up) {
    int gap = s / 2;

    while (gap > 0) {
        for (int i = 0; i + gap < s; ++i) {
            weave_swap(v, l + i, l + i + gap, up);
        }
        gap /= 2;
    }
}

void Sorter::iterative_weave(vector<int>& v, int n) {
    // Pad to next power of 2
    int power = 1;
    while (power < n) power *= 2;
    while (size(v) < power) v.push_back(INT_MAX);

    // Weave sort using bottom-up approach
    for (int s = 2; s <= power; s *= 2) {
        for (int i = 0; i < power; i += s) {
            // Inner block sort using weaving merge (bitonic-like pattern)
            iterative_weav(v, i, s, 1);
        }
    }

    // Remove padding
    while (size(v) > n) v.pop_back();
}

void Sorter::rotate_mer(vector<int>& v, int l, int m, int r) {
    // Find the point where the second part starts in sorted order
    int i = l;
    int j = m;

    while (i < j && j < r) {
        // If elements are already in order, skip
        if (v[i] <= v[j]) {
            ++i;
        }
        else {
            // Rotate the smaller element from right half into position
            int value = v[j];
            int index = j;

            // Shift elements from i to j-1 to the right by one
            while (index > i) {
                v[index] = v[index - 1];
                --index;
            }
            v[i] = value;

            // Update pointers
            ++i;
            ++j;
            ++m;
        }
    }
}

void Sorter::rotate_merg(vector<int>& v, int l, int r) {
    if (r - l <= 1) return;

    int m = l + (r - l) / 2;

    rotate_merg(v, l, m);
    rotate_merg(v, m, r);
    rotate_mer(v, l, m, r);
}

void Sorter::rotate_merge(vector<int>& v, int n) {
    rotate_merg(v, 0, n);
}

void Sorter::block_rotate(vector<int>& v, int l, int m, int r) {
    int temp = 0;
    int p = (l + m) / 2;
    for(int i=l;i<p;i++){
        temp = v[i];
        v[i] = v[m-1-(i-l)];
        v[m-1-(i-l)] = temp;
    }
    p = (m + r) / 2;
    for(int i=m;i<p;i++){
        temp = v[i];
        v[i] = v[r-1-(i-m)];
        v[r-1-(i-m)] = temp;
    }
    p = (l + r) / 2;
    for(int i=l;i<p;i++){
        temp = v[i];
        v[i] = v[r-1-(i-l)];
        v[r-1-(i-l)] = temp;
    }
}

int Sorter::block_search(const vector<int>& v, int l, int r, int key) {
    while (l < r) {
        int mid = l + (r - l) / 2;
        if (v[mid] < key)
            l = mid + 1;
        else
            r = mid;
    }
    return l;
}

void Sorter::block_merge(vector<int>& v, int l, int m, int r) {
    if (l >= m || m >= r) return;

    int m1 = l;
    int m2 = m;

    while (m1 < m && m2 < r) {
        if (v[m1] <= v[m2]) {
            m1++;
        }
        else {
            // Find how many elements in the right side are less than v[mid1] using binary search
            int index = block_search(v, m2, r, v[m1]);

            // Rotate the smaller block into place
            block_rotate(v, m1, m2, index);

            // Advance all pointers
            int shift = index - m2;
            m1 += shift + 1;
            m2 = index;
            m += shift; // move midpoint
        }
    }
}

void Sorter::bloc(vector<int>& v, int l, int r) {
    if (r - l <= 1) return;

    int m = (l + r) / 2;
    bloc(v, l, m);
    bloc(v, m, r);
    block_merge(v, l, m, r);
}

void Sorter::block(vector<int>& v, int n) {
    bloc(v, 0, n);
}

void Sorter::iterative_block_merge(vector<int>& v, int l, int m, int r) {
    if (l >= m || m >= r) return;

    int m1 = l;
    int m2 = m;

    while (m1 < m && m2 < r) {
        if (v[m1] <= v[m2]) {
            m1++;
        }
        else {
            // Find how many elements in the right side are less than v[mid1] using binary search
            int i = block_search(v, m2, r, v[m1]);

            // Rotate the smaller block into place
            block_rotate(v, m1, m2, i);

            // Advance all pointers
            int shift = i - m2;
            m1 += shift + 1;
            m2 = i;
            m += shift; // move midpoint
        }
    }
}

void Sorter::iterative_block(vector<int>& v, int n) {
    for (int s = 1; s < n; s *= 2) {
        for (int l = 0; l < n; l += 2 * s) {
            int m = min(l + s, n);
            int r = min(l + 2 * s, n);
            iterative_block_merge(v, l, m, r);
        }
    }
}

void Sorter::wiki(vector<int>& v, int n) {
    bool swapped = 1;
    while (swapped){
        swapped = 0;
        for (int i = 0; i < n - 1; i++) {
            // Swap elements if they're out of order
            if (v[i] > v[i + 1]) {
                swap(v[i], v[i + 1]);
                swapped = 1;
            }
        }
    }
}

void Sorter::grail_sub(vector<int>& v) {
    int n = size(v);
    if (n <= 1) return;

    int pivot = v[0];

    // Divide elements by relation to pivot
    vector<int> L, R;

    for (int i = 1; i < n; i++) {
        if (v[i] <= pivot) {
            L.push_back(v[i]);
        }
        else {
            R.push_back(v[i]);
        }
    }

    // Recurse on both subarrays
    grail_sub(L);
    grail_sub(R);

    // Combine the sorted groups back into the original array
    int j = 0;
    for(int i=0;i<size(L);i++) {
        v[j++] = L[i];
    }
    v[j++] = pivot;
    for(int i=0;i<size(R);i++) {
        v[j++] = R[i];
    }
}

void Sorter::grail(vector<int>& v, int n) {
    grail_sub(v);
}

void Sorter::stooge_sub(vector<int>& v, int l, int r) {
    int m;

    if (r - l + 1 > 2)
    {
        m = (r - l + 1) / 3;
        stooge_sub(v, l, r - m);
        stooge_sub(v, l + m, r);
        stooge_sub(v, l, r - m);
    }

    if (v[r] < v[l]) swap(v[r], v[l]);
}

void Sorter::stooge(vector<int>& v, int n)
{
    stooge_sub(v, 0, n-1);
}

void Sorter::weak_heap_sift(vector<int>& v, int i, int n) {
    int l = 2 * i + 1;
    int r = 2 * i + 2;
    int max = i;

    // Find the largest among root, l, r
    if (l < n && v[l] > v[max]) {
        max = l;
    }
    if (r < n && v[r] > v[max]) {
        max = r;
    }

    // If the max is not the root, swap and continue sifting down
    if (max != i) {
        swap(v[i], v[max]);
        weak_heap_sift(v, max, n);
    }
}

void Sorter::weak_heap_build(vector<int>& v, int n) {
    // Sift down elements starting from the last non-leaf node
    for (int i = n / 2 - 1; i >= 0; --i) {
        weak_heap_sift(v, i, n);
    }
}

void Sorter::weak_heap(vector<int>& v, int n) {
    weak_heap_build(v, n);

    // Extract elements from the heap one by one
    for (int i = n - 1; i > 0; --i) {
        // Swap the root (maximum element) with the last element
        swap(v[0], v[i]);

        // Restore the weak heap property
        weak_heap_sift(v, 0, i);
    }
}

int Sorter::leonardo(int k) {
    if (k < 2) return 1;
    return leonardo(k - 1) + leonardo(k - 2) + 1;
}

void Sorter::build_heap(vector<int>& v, int l, int r) {
    int i = l;
    int j = 0;
    int k = 0;

    while (k < r - l + 1) {
        if (k & 0xAAAAAAAA) {
            j = j + i;
            i = i >> 1;
        }
        else {
            i = i + j;
            j = j >> 1;
        }

        k = k + 1;
    }

    while (i > 0) {
        j = j >> 1;
        k = i + j;
        while (k < r) {
            if (v[k] > v[k - i]) {
                break;
            }
            swap(v[k], v[k - i]);
            k = k + i;
        }

        i = j;
    }
}

void Sorter::smooth(vector<int>& v, int n){

    int p = n - 1;
    int q = p;
    int r = 0;

    while (p > 0) {
        if ((r & 0x03) == 0) {
            build_heap(v, r, q);
        }

        if (leonardo(r) == p) {
            r = r + 1;
        }
        else {
            r = r - 1;
            q = q - leonardo(r);
            build_heap(v, r, q);
            q = r - 1;
            r = r + 1;
        }

        swap(v[0], v[p]);
        p--;
    }

    for (int i = 0; i < n - 1; i++) {
        int j = i + 1;
        while (j > 0 && v[j] < v[j - 1]) {
            swap(v[j], v[j - 1]);
            j = j - 1;
        }
    }
}

void Sorter::poplar_heap_build(vector<int>& v, int n, int i) {
    int max = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;

    // If left child is larger than root
    if (l < n && v[l] > v[max]) {
        max = l;
    }

    // If right child is larger than largest so far
    if (r < n && v[r] > v[max]) {
        max = r;
    }

    // If largest is not root, swap and continue heapifying
    if (max != i) {
        swap(v[i], v[max]);
        poplar_heap_build(v, n, max);
    }
}

void Sorter::poplar_heap(vector<int>& v, int n) {
    // Build the heap (rearrange the array)
    for (int i = n / 2 - 1; i >= 0; i--) {
        poplar_heap_build(v, n, i);
    }

    // One by one extract elements from the heap
    for (int i = n - 1; i >= 1; i--) {
        // Move current root to end
        swap(v[0], v[i]);

        // Heapify reduced heap
        poplar_heap_build(v, i, 0);
    }
}

int Sorter::binary_partition(vector<int>& v, int l, int r) {
    int pivot = v[r];  // Choose the last element as the pivot
    int i = l - 1;      // Index of smaller element

    // Move all elements smaller than pivot to the left
    for (int j = l; j < r; j++) {
        if (v[j] <= pivot) {
            i++;
            swap(v[i], v[j]);
        }
    }

    // Place pivot in its correct position
    swap(v[i + 1], v[r]);
    return i + 1;
}

void Sorter::binary_qui(vector<int>& v, int l, int r) {
    if (l < r) {
        // Partition the array
        int pi = binary_partition(v, l, r);

        // Recursively sort both subarrays
        binary_qui(v, l, pi - 1);
        binary_qui(v, pi + 1, r);
    }
}

void Sorter::binary_quick(vector<int>& v, int n) {
    binary_qui(v, 0, n - 1);
}

void Sorter::pancake_flip(vector<int>& v, int r) {
    int l = 0;
    while (l < r) {
        // Swap elements l and r
        int temp = v[l];
        v[l] = v[r];
        v[r] = temp;
        // Move l and r
        l++;
        r--;
    }
}

void Sorter::pancake(vector<int>& v, int n){

    for (int i = n - 1; i > 0; i--) {

        // Find the index of the maximum element in the unsorted portion of the vector
        int maxI = 0;
        int maximum = v[0];
        for (int j = 1; j <= i; j++) {
            if (v[j] > maximum) {
                maxI = j;
                maximum = v[j];
            }
        }
        // If the maximum element is already at the end of the unsorted portion, move on to the next iteration
        if (maxI == i) {
            continue;
        }

        // Flip the portion of the vector up to the maximum element
        pancake_flip(v, maxI);
        // Flip the entire unsorted portion of the vector to move the maximum element to the end of the unsorted portion
        pancake_flip(v, i);
    }
}

void Sorter::cocktail(vector<int>& v, int n){
    bool swapped = 1;
    int l = 0;
    int r = n - 1;

    while (swapped) {
        // reset the swapped flag on entering the loop, because it might be true from a previous iteration
        swapped = 0;

        // loop from left to right, like bubble sort
        for (int i = l; i < r; ++i) {
            if (v[i] > v[i + 1]) {
                swap(v[i], v[i + 1]);
                swapped = 1;
            }
        }

        // if nothing moved, then array is sorted.
        if (!swapped) break;

        // otherwise, reset the swapped flag so that it can be used in the next stage
        swapped = 0;

        // move the end point back by one, because item at the end is in its rightful spot
        --r;

        // from right to left, doing the
        // same comparison as in the previous stage
        for (int i = r - 1; i >= l; --i) {
            if (v[i] > v[i + 1]) {
                swap(v[i], v[i + 1]);
                swapped = true;
            }
        }

        // increase the starting point, because the last stage would have moved the next smallest number to its rightful spot
        ++l;
    }
}

void Sorter::odd_even(vector<int>& v, int n){
    bool sorted = 0;

    while (!sorted){
        sorted = 1;

        // Perform Bubble sort on odd indexed element
        for (int i = 1; i <= n - 2; i += 2){
            if (v[i] > v[i + 1]){
                swap(v[i], v[i + 1]);
                sorted = 0;
            }
        }

        // Perform Bubble sort on even indexed element
        for (int i = 0; i <= n - 2; i += 2){
            if (v[i] > v[i + 1]){
                swap(v[i], v[i + 1]);
                sorted = 0;
            }
        }
    }
}


bool Sorter::circle_pass(vector<int>& v, int l, int r) {
    bool swapped = 0;

    if (l == r) return 0;

    // Storing the upper and lower bounds of the list
    int lo = l, hi = r;

    // Perform swaps if needed
    while (lo < hi) {
        if (v[lo] > v[hi]) {
            swap(v[lo], v[hi]);
            swapped = 1;
        }
        lo++;
        hi--;
    }

    // Special case for odd-sized list
    if (lo == hi && v[lo] > v[hi + 1]) {
        swap(v[lo], v[hi + 1]);
        swapped = 1;
    }

    // Check sub-lists recursively
    int m = (l + r) / 2;

    bool c1 = circle_pass(v, l, m);
    bool c2 = circle_pass(v, m + 1, r);
    return swapped || c1 || c2;
}

void Sorter::circle(vector<int>& v, int n) {
    // Keep calling circle_pass while there is a swap operation
    while (circle_pass(v, 0, n - 1)) {}
}

void Sorter::merge_insert_insertion(vector<int>& v, int p, int q) {
    for (int i = p; i < q; i++) {
        int temp = v[i + 1];
        int j = i + 1;
        while (j > p && v[j - 1] > temp) {
            v[j] = v[j - 1];
            j--;
        }
        v[j] = temp;
    }
}

void Sorter::merge_insert_merge(vector<int>& v, int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;
    vector<int> L(n1, 0);
    vector<int> R(n2, 0);
    for (int i = 0; i < n1; i++) {
        L.push_back(v[i + l]);
    }
    for (int i = 0; i < n2; i++) {
        R.push_back(v[i + m + 1]);
    }
    
    v.clear();

    int ri = 0;
    int li = 0;
    for (int i = l; i < r - l + 1; i++) {
        if (ri == n2) {
            v.push_back(L[li]);
            li++;
        }
        else if (li == n1) {
            v.push_back(R[ri]);
            ri++;
        }
        else if (R[ri] > L[li]) {
            v.push_back(L[li]);
            li++;
        }
        else {
            v.push_back(R[ri]);
            ri++;
        }
    }
}

void Sorter::merge_insert(vector<int>& v, int l, int r, int b) {
    if (r - l > b) {
        int m = (l + r) / 2;
        merge_insert(v, l, m, b);
        merge_insert(v, m + 1, r, b);
        merge_insert_merge(v, l, m, r);
    }
    else {
        merge_insert_insertion(v, l, r);
    }
}

void Sorter::merge_insertion(vector<int>& v, int n, int b) {
    merge_insert(v, 0, n-1, b);
}

struct Node {
    int v;
    Node* l;
    Node* r;

    Node(int v) : v(v), l(nullptr), r(nullptr) {}
};

Node* Sorter::tree_insert(Node* root, int v) {
    if (!root) return new Node(v);

    if (v <= root->v)
        root->l = tree_insert(root->l, v);
    else
        root->r = tree_insert(root->r, v);

    return root;
}

void Sorter::tree_traverse(Node* root, vector<int>& sorted) {
    if (!root) return;
    tree_traverse(root->l, sorted);
    sorted.push_back(root->v);
    tree_traverse(root->r, sorted);
}

void Sorter::tree(vector<int>& v, int n) {
    Node* root = nullptr;

    // Build binary search tree
    for (int i = 0; i < n; ++i) {
        root = tree_insert(root, v[i]);
    }

    // Get sorted values
    vector<int> sorted;
    tree_traverse(root, sorted);

    // Copy back to original vector
    v = sorted;
}

void Sorter::tournament_build(vector<int>& tree, vector<int>& v, int n) {
    for (int i = 0; i < n; ++i)
        tree[n + i] = v[i];

    for (int i = n - 1; i > 0; --i)
        tree[i] = min(tree[2 * i], tree[2 * i + 1]);
}

void Sorter::tournament_update(vector<int>& tree, int pos, int value, int n) {
    pos += n;
    tree[pos] = value;

    while (pos > 1) {
        pos /= 2;
        tree[pos] = min(tree[2 * pos], tree[2 * pos + 1]);
    }
}

int Sorter::tournament_min(vector<int>& tree, int n) {
    int min = tree[1];
    for (int i = 0; i < n; ++i) {
        if (tree[n + i] == min) return i;
    }
    return -1;
}

void Sorter::tournament(vector<int>& v, int n) {
    // Tree size must be 2n
    vector<int> tree(2 * n);

    // Build initial tree
    tournament_build(tree, v, n);

    // Repeatedly extract the min and update tree
    vector<int> sorted;
    for (int i = 0; i < n; ++i) {
        int minI = tournament_min(tree, n);
        sorted.push_back(v[minI]);
        tournament_update(tree, minI, INT_MAX, n); // Replace with INT_MAX
    }

    // Copy sorted result back to v
    v = sorted;
}

typedef pair<int, int> Pair;

void Sorter::stable_tournament_build(vector<Pair>& tree, const vector<Pair>& v, int n) {
    for (int i = 0; i < n; ++i)
        tree[n + i] = v[i];

    for (int i = n - 1; i > 0; --i)
        tree[i] = min(tree[2 * i], tree[2 * i + 1]);
}

void Sorter::stable_tournament_update(vector<Pair>& tree, int pos, Pair value, int n) {
    pos += n;
    tree[pos] = value;

    while (pos > 1) {
        pos /= 2;
        tree[pos] = min(tree[2 * pos], tree[2 * pos + 1]);
    }
}

int Sorter::stable_tournament_min(const vector<Pair>& tree, int n, Pair minVal) {
    for (int i = 0; i < n; ++i) {
        if (tree[n + i] == minVal)
            return i;
    }
    return -1;
}

void Sorter::stable_tournament(vector<int>& v, int n) {
    // Store original indices to maintain stability
    vector<Pair> original(n);
    for (int i = 0; i < n; ++i)
        original[i] = { v[i], i };

    // Tree size = 2n
    vector<Pair> tree(2 * n, { INT_MAX, INT_MAX });

    stable_tournament_build(tree, original, n);

    vector<int> sorted;
    for (int i = 0; i < n; ++i) {
        Pair min = tree[1];
        int minI = stable_tournament_min(tree, n, min);
        sorted.push_back(min.first);
        stable_tournament_update(tree, minI, { INT_MAX, INT_MAX }, n);
    }

    // Copy sorted result back to v
    v = sorted;
}

void Sorter::gnome(vector<int>& v, int n){
    int i = 0;

    while (i < n) {
        if (i == 0)
            i++;
        if (v[i] >= v[i - 1])
            i++;
        else {
            swap(v[i], v[i - 1]);
            i--;
        }
    }
}

void Sorter::library(vector<int>& v, int n) {
    if (n <= 1) return;

    double numGaps = 1.0;

    // Initial size of the table with gaps
    int size = (int)(ceil((1.0 + numGaps) * (double)n));
    vector<int> table(size, INT_MIN);

    // Insert the first element
    int count = 1;
    table[size / 2] = v[0];

    for (int i = 1; i < n; i++) {
        int l = 0, r = size - 1;
        int key = v[i];

        // Binary search to find the right spot
        while (l <= r) {
            int m = (l + r) / 2;

            // Skip empty spots
            while (m <= r && table[m] == INT_MIN) m++;
            if (m > r || table[m] > key) {
                r = min(m - 1, r - 1);
            }else {
                l = m + 1;
            }
        }

        // Shift elements if needed to make space
        int pos = l;
        while (pos < size && table[pos] != INT_MIN) pos++;
        
        if (pos >= size) {
            // Resize and reinsert if out of space
            size = (int)(ceil((1.0 + numGaps) * (i + 1)));
            vector<int> newTable(size, INT_MIN);

            int j = 0;
            for (int k = 0; k < size && j < count; k++) {
                if (k % (int)(numGaps + 1.0) == 0) {
                    newTable[k] = table[j];
                    j++;
                }
            }

            table = newTable;

            // Redo insertion
            i--;
            continue;
        }

        table[pos] = key;
        count++;
        //for(int j=0;j<table.size();j++){
        //    cout << table[j] << " ";
        //}
        //cout << "\n";
    }

    // Extract sorted elements back to the original vector
    v.clear();
    for (int i = 0; i < size; i++) {
        if (table[i] != INT_MIN) v.push_back(table[i]);
    }
}

vector<int> Sorter::strand_merge(vector<int>& a, vector<int>& b) {
    vector<int> out;
    auto itA = a.begin();
    auto itB = b.begin();

    while (itA != a.end() && itB != b.end()) {
        if (*itA < *itB) {
            out.push_back(*itA++);
        }
        else {
            out.push_back(*itB++);
        }
    }

    while (itA != a.end()) out.push_back(*itA++);
    while (itB != b.end()) out.push_back(*itB++);

    return out;
}

void Sorter::strand(vector<int>& v, int n) {
    vector<int> in(v.begin(), v.end());
    vector<int> out;

    while (!in.empty()) {
        vector<int> strand;
        auto it = in.begin();
        strand.push_back(*it);
        it = in.erase(it);

        while(it != in.end()) {
            if (*it >= strand.back()) {
                strand.push_back(*it);
                it = in.erase(it);
            }
            else {
                ++it;
            }
        }

        out = strand_merge(out, strand);
    }

    // Copy sorted result back to the vector
    v.assign(out.begin(), out.end());
}

void Sorter::patience_merge(vector<int>& v, vector<vector<int>>& p)
{
    vector<int> out;

    // Always move smallest (top) to out
    while (1) {

        // Stores the smallest element of the top of the piles
        int min = INT_MAX;

        // Stores index of the smallest element of the top of the piles
        int index = -1;

        // Calculate the smallest element
        for (int i = 0; i < size(p); i++) {

            // If min is greater than top of the current stack
            if (min > p[i][size(p[i]) - 1]) {

                // Update min
                min = p[i][size(p[i]) - 1];

                // Update index
                index = i;
            }
        }

        // Insert the smallest element of the top of the stack
        out.push_back(min);

        // Remove the top element from the current pile
        p[index].pop_back();

        if (p[index].empty()) {

            // Remove current pile from all piles
            p.erase(p.begin() + index);
        }

        // If all the piles are empty
        if (p.size() == 0) break;
    }

    v = out;
}

void Sorter::patience(vector<int>& v, int n){
    vector<vector<int>> p;

    for (int i = 0; i < n; i++) {

        // If no piles are created
        if (p.empty()) {

            // Initialize a new pile in p with current element
            vector<int> temp = { v[i] };
            p.push_back(temp);
        }
        else {

            // Check if top element of each pile is less than or equal to current element or not
            int flag = 1;

            // Traverse piles
            for (int j = 0; j < size(p); j++) {

                // Check if the element to be inserted is less than current pile's top
                if (v[i] < p[j][size(p[j]) - 1]) {
                    p[j].push_back(v[i]);

                    // Update flag
                    flag = 0;
                    break;
                }
            }

            // If flag is true
            if (flag) {

                // Initialize a new pile in p with current element
                vector<int> temp = { v[i] };
                p.push_back(temp);
            }
        }
    }

    // Sort the array
    patience_merge(v, p);
}

void Sorter::bitonic_merge(vector<int>& v, int l, int c, bool dir) {
    if (c <= 1) return;

    int k = c / 2;
    for (int i = l; i < l + k; i++){
        if (dir == v[i] > v[i + k]){
            swap(v[i], v[i + k]);
        }
    }

    bitonic_merge(v, l, k, dir);
    bitonic_merge(v, l + k, k, dir);
}

void Sorter::bitonic_sort(vector<int>& v, int l, int c, bool dir) {
    if (c <= 1) return;

    int k = c / 2;

    // Sort first half ascending and second half descending
    bitonic_sort(v, l, k, 1);
    bitonic_sort(v, l + k, k, 0);

    // Merge whole sequence in desired direction
    bitonic_merge(v, l, c, dir);
}

void Sorter::bitonic(vector<int>& v, int n) {
    
    int power = 1;
    while (power < n) power *= 2;
    for (int i = n; i < power; i++) v.push_back(INT_MAX); // Padding

    bitonic_sort(v, 0, power, 1);

    // Remove padding
    v.resize(n);
}

void Sorter::odd_even_mer(vector<int>& v, int l, int n, int r) {
    int step = r * 2;
    if (step < n) {
        odd_even_mer(v, l, n, step);
        odd_even_mer(v, l + r, n, step);
        for (int i = l + r; i + r < l + n; i += step) {
            if (v[i] > v[i + r]) {
                swap(v[i], v[i + r]);
            }
        }
    } else {
        if (l + r < l + n && v[l] > v[l + r]) {
            swap(v[l], v[l + r]);
        }
    }
}

void Sorter::odd_even_sort(vector<int>& v, int l, int n) {
    if (n > 1) {
        int m = n / 2;
        odd_even_sort(v, l, m);
        odd_even_sort(v, l + m, m);
        odd_even_mer(v, l, n, 1);
    }
}

void Sorter::odd_even_merge(vector<int>& v, int n) {
    // Ensure n is a power of two
    int power = 1;
    while (power < n) power *= 2;

    // Pad with INT_MAX if needed
    for (int i = n; i < power; i++)
        v.push_back(INT_MAX);

    odd_even_sort(v, 0, power);

    // Remove padding
    v.resize(n);
}

void Sorter::pairwise_network_merge(vector<int>& v, int l, int r) {
    int m = (l + r) / 2;
    if (r - l <= 1) return;
    int h = m - l;

    for (int i = l; i < m; i++) {
        if (v[i] > v[i + h]) {
            swap(v[i], v[i + h]);
        }
    }

    pairwise_network_merge(v, l, m);
    pairwise_network_merge(v, m, r);
}

void Sorter::pairwise_network_sort(vector<int>& v, int l, int r) {
    if (r - l <= 1) return;

    int m = (l + r) / 2;
    pairwise_network_sort(v, l, m);
    pairwise_network_sort(v, m, r);
    pairwise_network_merge(v, l, r);
}

void Sorter::pairwise_network(vector<int>& v, int n) {
    int a = 1;
    while(a < n){
        int b = a;
        int c = 0;
        while(b < n){
            if(v[b-a] > v[b]) swap(v[b-a], v[b]);
            b++;
            c++;
            if(c >= a){
                c = 0;
                b += a;
            }
        }
    a *= 2;
    }
    a /= 4;
    int e = 1;
    while(a > 0){
        int d = e;
        while(d > 0){
            int b = (d + 1) * a;
            int c = 0;
            while(b < n){
                if(v[b-d*a] > v[b]) swap(v[b-d*a], v[b]);
                b++;
                c++;
                if(c >= a){
                    c = 0;
                    b += a;
                }
            }
            d /= 2;
        }
        a /= 2;
        e *= 2;
        e++;
    }
}

int Sorter::spread_digits(int value) {
    int pos = 0;
    while (value >>= 1) ++pos;
    return pos;
}

void Sorter::spread_sort(vector<int>& v, int l, int r, int shift, int b) {
    if (l >= r || shift < 0) return;

    const int bucketCount = 1 << b;
    vector<vector<int>> buckets(bucketCount);

    // Partition elements into buckets based on current bits
    for (int i = l; i <= r; ++i) {
        int bucketIndex = (v[i] >> shift) & (bucketCount - 1);
        buckets[bucketIndex].push_back(v[i]);
    }

    // Recursively sort each bucket and concatenate
    int index = l;
    int temp = 0;
    for (int i = 0; i < bucketCount; ++i) {
        if (!buckets[i].empty()) {
            if (size(buckets[i]) > bucketCount) {
                // Recursively sort large buckets
                for(int val: buckets[i]) v[index++] = val;
                spread_sort(v, index - size(buckets[i]), index - 1, shift - b, b);
            }
            else {
                // Sort for small buckets
                int s = size(buckets[i]);
                for(int j=0;j<s;j++){
                    for(int k=j+1;k<s;k++){
                        if(buckets[i][j] > buckets[i][k]){
                            temp = buckets[i][j];
                            buckets[i][j] = buckets[i][k];
                            buckets[i][k] = temp;
                        }
                    }
                }
                for(int val: buckets[i]) v[index++] = val;
            }
        }
    }
}

void Sorter::spread(vector<int>& v, int n, int b) {
    if (n == 0) return;

    int maxVal = v[0];
    for(int i=1;i<n;i++){
        maxVal = max(maxVal, v[i]);
    }
    int shift = spread_digits(maxVal);

    // Round up to nearest multiple of b bits
    shift = (shift / b) * b;
    spread_sort(v, 0, n - 1, shift, b);
}

void Sorter::flash(vector<int>& v, int n, float b) {
    if (n <= 1) return;

    int min = v[0], maxI = 0;
    for (int i = 1; i < n; ++i) {
        if (v[i] < min) min = v[i];
        if (v[i] > v[maxI]) maxI = i;
    }

    if (min == v[maxI]) return; // Already sorted (all elements equal)

    const int m = int(b * n);  // Customizable number of classes, typically 0.43
    vector<int> l(m, 0);

    // Classification
    double c = double(m - 1) / (v[maxI] - min);
    for (int i = 0; i < n; ++i)
        ++l[int(c * (v[i] - min))];

    // Accumulate counts
    for (int i = 1; i < m; ++i)
        l[i] += l[i - 1];

    // Permutation
    swap(v[0], v[maxI]); // Move max to front
    int move = 0, j = 0, k = m - 1;
    while (move < n - 1) {
        while (j > l[k] - 1)
            k = int(c * (v[++j] - min));

        int remove = v[j];
        while (j != l[k]) {
            k = int(c * (remove - min));
            int temp = v[l[k] - 1];
            v[l[k] - 1] = remove;
            remove = temp;
            --l[k];
            ++move;
        }
    }

    // Final cleanup: insertion sort
    for (int i = 1; i < n; ++i) {
        int key = v[i];
        int j = i - 1;
        while (j >= 0 && v[j] > key) {
            v[j + 1] = v[j];
            --j;
        }
        v[j + 1] = key;
    }
}

void Sorter::bead(vector<int>& v, int n) {
    // Find the maximum element 
    int maxi = v[0];
    for (int i = 1; i < n; i++) {
        if (v[i] > maxi) {
            maxi = v[i];
        }
    }

    // create empty memory 
    vector<vector<int>> beads;
    beads.resize(n);
    for (int i = 0; i < n; i++) {
        beads[i].resize(maxi);
        fill(beads[i].begin(), beads[i].end(), 0);
    }

    // mark the beads 
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < v[i]; j++) {
            beads[i][j] = 1;
        }
    }

    // move down the beads 
    for (int j = 0; j < maxi; j++) {
        int sum = 0;
        for (int i = 0; i < n; i++) {
            sum += beads[i][j];
            beads[i][j] = 0;
        }
        for (int i = n - 1; i >= n - sum; i--) {
            beads[i][j] = 1;
        }
    }

    // to get sorted array back into v
    for (int i = 0; i < n; i++) {
        int sum = 0;
        for (int j = 0; j < maxi; j++) {
            sum += beads[i][j];
        }
        v[i] = sum;
    }
}

void Sorter::spaghetti(vector<int>& v, int n) {

    int maxVal = v[0];
    for (int i=1;i<n;i++) {
        maxVal = max(maxVal, v[i]);
    }

    vector<int> rods(maxVal + 1, 0);

    for(int i=0;i<n;i++) rods[v[i]]++;

    // Collect spaghetti in order
    int j = 0;
    for(int i=0;i<=maxVal;i++) {
        while (rods[i]--) {
            v[j++] = i;
        }
    }
}

void Sorter::ska_digit(vector<int>& v, vector<int>& temp, int byte) {
    int count[256] = { 0 };

    // Count byte occurrences
    for (int val: v) {
        int b = (val >> (byte * 8)) & 0xFF;
        ++count[b];
    }

    // Compute bucket starts
    int offset[256];
    offset[0] = 0;
    for (int i = 1; i < 256; ++i) offset[i] = offset[i - 1] + count[i - 1];

    // Sort based on byte
    for (int val: v) {
        int b = (val >> (byte * 8)) & 0xFF;
        temp[offset[b]++] = val;
    }

    // Copy back
    v.swap(temp);
}

void Sorter::ska(vector<int>& v, int n) {
    if (n <= 1) return;

    vector<int> temp(n);

    for (int byte = 0; byte < 4; ++byte) {
        ska_digit(v, temp, byte);
    }
}

void fillVector() {
    v.clear();

    for (int i = 0; i < n; i++) {
        v.push_back(rand());
    }
}

bool checkForSorted() {
    for (int i = 1; i < n; i++) {
        if (v[i] < v[i - 1]) return 0;
    }
    return 1;
}

void sortVector(int index, Sorter sorter, vector<int>& v, int n) {
    switch (index) {
    case 0: // bubble
        sorter.bubble(v, n);
        break;
    case 1: // exchange
        sorter.exchange(v, n);
        break;
    case 2: // selection
        sorter.selection(v, n);
        break;
    case 3: // double selection
        sorter.double_selection(v, n);
        break;
    case 4: // insertion
        sorter.insertion(v, n);
        break;
    case 5: // binary insertion
        sorter.binary_insertion(v, n);
        break;
    case 6: // merge
        sorter.merge(v, n);
        break;
    case 7: // three-way merge
        sorter.merge_3(v, n);
        break;
    case 8: // four-way merge
        sorter.merge_4(v, n);
        break;
    case 9: // quick, lomuto partition
        sorter.quick_lomuto(v, n);
        break;
    case 10: // quick, hoare partition
        sorter.quick_hoare(v, n);
        break;
    case 11: // quick, naive partition
        sorter.quick_naive(v, n);
        break;
    case 12: // heap
        sorter.heap(v, n);
        break;
    case 13: // comb
        sorter.comb(v, n);
        break;
    case 14: // shell
        sorter.shell(v, n);
        break;
    case 15: // bogo
        sorter.bogo(v, n);
        break;

    case 16: // pigeonhole
        sorter.pigeonhole(v, n);
        break;
    case 17: // counting
        sorter.counting(v, n);
        break;
    case 18: // bucket, 4 buckets
        sorter.bucket(v, n, 4);
        break;
    case 19: // bucket, 16 buckets
        sorter.bucket(v, n, 16);
        break;
    case 20: // bucket, 64 buckets
        sorter.bucket(v, n, 64);
        break;
    case 21: // bucket, N buckets
        sorter.bucket(v, n, n);
        break;
    case 22: // lsd radix, base 2
        sorter.lsd(v, n, 2);
        break;
    case 23: // lsd radix, base 4
        sorter.lsd(v, n, 4);
        break;
    case 24: // lsd radix, base 8
        sorter.lsd(v, n, 8);
        break;
    case 25: // lsd radix, base 16
        sorter.lsd(v, n, 16);
        break;
    case 26: // msd radix, base 2
        sorter.msd(v, n, 2);
        break;
    case 27: // msd radix, base 4
        sorter.msd(v, n, 4);
        break;
    case 28: // msd radix, base 8
        sorter.msd(v, n, 8);
        break;
    case 29: // msd radix, base 16
        sorter.msd(v, n, 16);
        break;
    case 30: // in place lsd radix, base 2
        sorter.in_place_lsd(v, n, 2);
        break;
    case 31: // in place lsd radix, base 4
        sorter.in_place_lsd(v, n, 4);
        break;
    case 32: // in place lsd radix, base 8
        sorter.in_place_lsd(v, n, 8);
        break;
    case 33: // in place lsd radix, base 16
        sorter.in_place_lsd(v, n, 16);
        break;
    case 34: // in place msd radix, base 2
        sorter.in_place_msd(v, n, 2);
        break;
    case 35: // in place msd radix, base 4
        sorter.in_place_msd(v, n, 4);
        break;
    case 36: // in place msd radix, base 8
        sorter.in_place_msd(v, n, 8);
        break;
    case 37: // in place msd radix, base 16
        sorter.in_place_msd(v, n, 16);
        break;

    case 38: // intro
        sorter.intro(v, n);
        break;
    case 39: // cycle
        sorter.cycle(v, n);
        break;
    case 40: // tim, group size 4
        sorter.tim(v, n, 4);
        break;
    case 41: // tim, group size 16
        sorter.tim(v, n, 16);
        break;
    case 42: // tim, group size 64
        sorter.tim(v, n, 64);
        break;
    case 43: // iterative merge
        sorter.iterative_merge(v, n);
        break;
    case 44: // naive in-place merge
        sorter.in_place_merge(v, n);
        break;
    case 45: // weave
        sorter.weave(v, n);
        break;
    case 46: // iterative weave
        sorter.iterative_weave(v, n);
        break;
    case 47: // rotate merge
        sorter.rotate_merge(v, n);
        break;
    case 48: // block
        sorter.block(v, n);
        break;
    case 49: // iterative block
        sorter.iterative_block(v, n);
        break;
    case 50: // wiki
        sorter.wiki(v, n);
        break;
    case 51: // grail
        sorter.grail(v, n);
        break;
    case 52: // stooge
        sorter.stooge(v, n);
        break;
    case 53: // weak heap
        sorter.weak_heap(v, n);
        break;
    case 54: // smooth
        sorter.smooth(v, n);
        break;
    case 55: // poplar heap
        sorter.poplar_heap(v, n);
        break;
    case 56: // binary quick
        sorter.binary_quick(v, n);
        break;
    case 57: // pancake
        sorter.pancake(v, n);
        break;
    case 58: // cocktail shaker
        sorter.cocktail(v, n);
        break;
    case 59: // odd-even
        sorter.odd_even(v, n);
        break;
    case 60: // circle
        sorter.circle(v, n);
        break;
    case 61: // merge-insertion size 8
        sorter.merge_insertion(v, n, 8);
        break;
    case 62: // tree
        sorter.tree(v, n);
        break;
    case 63: // tournament
        sorter.tournament(v, n);
        break;
    case 64: // stable tournament
        sorter.stable_tournament(v, n);
        break;
    case 65: // gnome
        sorter.gnome(v, n);
        break;
    case 66: // library
        sorter.library(v, n);
        break;
    case 67: // strand
        sorter.strand(v, n);
        break;
    case 68: // patience
        sorter.patience(v, n);
        break;
    case 69: // bitonic
        sorter.bitonic(v, n);
        break;
    case 70: // odd-even merge
        sorter.odd_even_merge(v, n);
        break;
    case 71: // pairwise network
        sorter.pairwise_network(v, n);
        break;

    case 72: // spread, 1 bit base
        sorter.spread(v, n, 1);
        break;
    case 73: // spread, 2 bit base
        sorter.spread(v, n, 2);
        break;
    case 74: // spread, 4 bit base
        sorter.spread(v, n, 4);
        break;
    case 75: // spread, 8 bit base
        sorter.spread(v, n, 8);
        break;
    case 76: // flash, 0.3n buckets
        sorter.flash(v, n, 0.3);
        break;
    case 77: // flash, 0.4n buckets
        sorter.flash(v, n, 0.4);
        break;
    case 78: // flash, 0.5n buckets
        sorter.flash(v, n, 0.5);
        break;
    case 79: // bead (gravity)
        sorter.bead(v, n);
        break;
    case 80: // spaghetti
        sorter.spaghetti(v, n);
        break;
    case 81: // ska
        sorter.ska(v, n);
        break;
    }
}

void runAllSorts(bool print) {

    Sorter sorter = Sorter();

    for(int i = 0; i < numSorts; i++){
        // skip bogo
        if(i == 15) continue;
        // skip library
        if(i == 66) continue;
        
        fillVector();

        auto start = chrono::high_resolution_clock::now();

        sorter.sortVector(i, v, n);

        auto end = chrono::high_resolution_clock::now();

        chrono::nanoseconds diff = end - start;
        int d = (int)diff.count();

        if(print){
            if(checkForSorted()){
                cout << "Sort #" << i + 1 << " (" << sorts[i].name << ") succeeded in " << d << " ns.\n";
            }else{
                cout << "Sort #" << i + 1 << " (" << sorts[i].name << ") failed in " << d << " ns.\n";
            }
        }
    }
}

void initAllSorts() {

    // POPULAR COMPARISON SORTS
    sorts[0] = Sort(0, "Bubble", 1, "Exchange", 1, 1, 1, "n", "n^2", "n^2", "c");
    sorts[1] = Sort(1, "Exchange", 1, "Exchange", 0, 1, 1, "n^2", "n^2", "n^2", "c");
    sorts[2] = Sort(2, "Selection", 1, "Selection", 0, 1, 1, "n^2", "n^2", "n^2", "c");
    sorts[3] = Sort(3, "Double Selection", 1, "Selection", 0, 1, 1, "n^2", "n^2", "n^2", "c");
    sorts[4] = Sort(4, "Insertion", 1, "Insertion", 1, 1, 1, "n", "n^2", "n^2", "c");
    sorts[5] = Sort(5, "Binary Insertion", 1, "Insertion", 1, 1, 1, "n log n", "n^2", "n^2", "c");
    sorts[6] = Sort(6, "Merge", 1, "Merge", 1, 0, 1, "n log n", "n log n", "n log n", "n");
    sorts[7] = Sort(7, "Three-Way Merge", 1, "Merge", 1, 0, 1, "n log n", "n log n", "n log n", "n");
    sorts[8] = Sort(8, "Four-Way Merge", 1, "Merge", 1, 0, 1, "n log n", "n log n", "n log n", "n");
    sorts[9] = Sort(9, "Quick, Lomuto Partition", 1, "Partitioning", 0, 0, 1, "n log n", "n log n", "n^2", "log n");
    sorts[10] = Sort(10, "Quick, Hoare Partition", 1, "Partitioning", 0, 0, 1, "n log n", "n log n", "n^2", "log n");
    sorts[11] = Sort(11, "Quick, Naive Partition", 1, "Partitioning", 0, 0, 1, "n log n", "n log n", "n^2", "log n");
    sorts[12] = Sort(12, "Heap", 1, "Selection", 0, 1, 1, "n log n", "n log n", "n log n", "c");
    sorts[13] = Sort(13, "Comb", 1, "Exchange", 0, 1, 1, "n log n", "n^1.2", "n^2", "c");
    sorts[14] = Sort(14, "Shell", 1, "Insertion", 0, 1, 1, "n log n", "n^1.2", "n^1.5", "c");
    sorts[15] = Sort(15, "Bogo", 1, "Exchange", 0, 1, 1, "n", "n!n", "infinite", "c");

    // POPULAR NON-COMPARISON SORTS
    sorts[16] = Sort(16, "Pigeonhole", 0, "Buckets", 1, 0, 0, "", "", "", "");
    sorts[17] = Sort(17, "Counting", 0, "Buckets", 1, 0, 0, "n+b^d", "n+b^d", "n+b^d", "n+b^d");
    sorts[18] = Sort(18, "Bucket, 4 Buckets", 0, "Buckets", 1, 0, 1, "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b");
    sorts[19] = Sort(19, "Bucket, 16 Buckets", 0, "Buckets", 1, 0, 1, "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b");
    sorts[20] = Sort(20, "Bucket, 64 Buckets", 0, "Buckets", 1, 0, 1, "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b");
    sorts[21] = Sort(21, "Bucket, N Buckets", 0, "Buckets", 1, 0, 1, "3n", "3n", "3n", "2n");
    sorts[22] = Sort(22, "LSD Radix, Base 2", 0, "Buckets", 1, 0, 1, "nd+bd", "nd+bd", "nd+bd", "n+b");
    sorts[23] = Sort(23, "LSD Radix, Base 4", 0, "Buckets", 1, 0, 1, "nd+bd", "nd+bd", "nd+bd", "n+b");
    sorts[24] = Sort(24, "LSD Radix, Base 8", 0, "Buckets", 1, 0, 1, "nd+bd", "nd+bd", "nd+bd", "n+b");
    sorts[25] = Sort(25, "LSD Radix, Base 16", 0, "Buckets", 1, 0, 1, "nd+bd", "nd+bd", "nd+bd", "n+b");
    sorts[26] = Sort(26, "MSD Radix, Base 2", 0, "Buckets", 1, 0, 1, "ndb", "ndb", "ndb", "n+b+d");
    sorts[27] = Sort(27, "MSD Radix, Base 4", 0, "Buckets", 1, 0, 1, "ndb", "ndb", "ndb", "n+b+d");
    sorts[28] = Sort(28, "MSD Radix, Base 8", 0, "Buckets", 1, 0, 1, "ndb", "ndb", "ndb", "n+b+d");
    sorts[29] = Sort(29, "MSD Radix, Base 16", 0, "Buckets", 1, 0, 1, "ndb", "ndb", "ndb", "n+b+d");
    sorts[30] = Sort(30, "In-Place LSD Radix, Base 2", 0, "Buckets", 1, 1, 1, "dn^2", "dn^2", "dn^2", "b");
    sorts[31] = Sort(31, "In-Place LSD Radix, Base 4", 0, "Buckets", 1, 1, 1, "dn^2", "dn^2", "dn^2", "b");
    sorts[32] = Sort(32, "In-Place LSD Radix, Base 8", 0, "Buckets", 1, 1, 1, "dn^2", "dn^2", "dn^2", "b");
    sorts[33] = Sort(33, "In-Place LSD Radix, Base 16", 0, "Buckets", 1, 1, 1, "dn^2", "dn^2", "dn^2", "b");
    sorts[34] = Sort(34, "In-Place MSD Radix, Base 2", 0, "Buckets", 1, 1, 1, "ndb", "ndb", "ndb", "b+d");
    sorts[35] = Sort(35, "In-Place MSD Radix, Base 4", 0, "Buckets", 1, 1, 1, "ndb", "ndb", "ndb", "b+d");
    sorts[36] = Sort(36, "In-Place MSD Radix, Base 8", 0, "Buckets", 1, 1, 1, "ndb", "ndb", "ndb", "b+d");
    sorts[37] = Sort(37, "In-Place MSD Radix, Base 16", 0, "Buckets", 1, 1, 1, "ndb", "ndb", "ndb", "b+d");

    // OBSCURE COMPARISON SORTS
    sorts[38] = Sort(38, "Intro", 1, "Partitioning, Selection", 0, 0, 1, "n log n", "n log n", "n log n", "log n");
    sorts[39] = Sort(39, "Cycle", 1, "Selection", 0, 1, 1, "n^2", "n^2", "n^2", "c");
    sorts[40] = Sort(40, "Tim, 4", 1, "Merge", 1, 0, 1, "n", "n log n", "n log n", "n");
    sorts[41] = Sort(41, "Tim, 16", 1, "Merge", 1, 0, 1, "n", "n log n", "n log n", "n");
    sorts[42] = Sort(42, "Tim, 64", 1, "Merge", 1, 0, 1, "n", "n log n", "n log n", "n");
    sorts[43] = Sort(43, "Iterative Merge", 1, "Merge", 1, 0, 1, "n log n", "n log n", "n log n", "n");
    sorts[44] = Sort(44, "Naive In-Place Merge", 1, "Merge", 1, 1, 1, "n^2", "n^2", "n^2", "c");
    sorts[45] = Sort(45, "Weave", 1, "Merge", 0, 1, 1, "n^2", "n^2", "n^2", "c");
    sorts[46] = Sort(46, "Iterative Weave", 1, "Merge", 0, 1, 1, "n^2", "n^2", "n^2", "c");
    sorts[47] = Sort(47, "Rotate Merge", 1, "Merge", 1, 1, 1, "n log n", "n^2", "n^2", "c");
    sorts[48] = Sort(48, "Block", 1, "Merge", 1, 1, 1, "n log n", "n log n", "n log n", "c");
    sorts[49] = Sort(49, "Iterative Block", 1, "Merge", 1, 1, 1, "n log n", "n log n", "n log n", "c");
    sorts[50] = Sort(50, "Wiki", 1, "Exchange", 1, 1, 1, "n", "n^2", "n^2", "c");
    sorts[51] = Sort(51, "Grail", 1, "Merge", 0, 0, 1, "n log n", "n log n", "n^2", "n");
    sorts[52] = Sort(52, "Stooge", 1, "Exchange", 0, 0, 1, "n", "n^2.7", "n^2.7", "n");
    sorts[53] = Sort(53, "Weak Heap", 1, "Exchange", 0, 0, 1, "n log n", "n log n", "n log n", "n");
    sorts[54] = Sort(54, "Smooth", 1, "Exchange", 0, 1, 1, "n", "n log n", "n log n", "c");
    sorts[55] = Sort(55, "Poplar Heap", 1, "Exchange", 0, 1, 1, "n log n", "n log n", "n log n", "c");
    sorts[56] = Sort(56, "Binary Quick", 1, "Exchange", 0, 0, 1, "n log n", "n log n", "n^2", "log n");
    sorts[57] = Sort(57, "Pancake", 1, "Exchange", 0, 1, 1, "n", "n^2", "n^2", "c");
    sorts[58] = Sort(58, "Cocktail Shaker", 1, "Exchange", 1, 1, 1, "n", "n^2", "n^2", "c");
    sorts[59] = Sort(59, "Odd-Even", 1, "Exchange", 1, 1, 1, "n", "n^2", "n^2", "c");
    sorts[60] = Sort(60, "Circle", 1, "Exchange", 1, 0, 1, "n log n", "n log n log n", "n log n log n", "log n");
    sorts[61] = Sort(61, "Merge-Insertion", 1, "Merge, Insertion", 0, 0, 1, "n", "n^2", "n^2", "log n");
    sorts[62] = Sort(62, "Tree", 1, "Insertion", 1, 0, 1, "n log n", "n log n", "n log n", "n");
    sorts[63] = Sort(63, "Tournament", 1, "Selection", 0, 0, 1, "n log n", "n log n", "n log n", "n");
    sorts[64] = Sort(64, "Stable Tournament", 1, "Selection", 1, 0, 1, "n log n", "n log n", "n log n", "n");
    sorts[65] = Sort(65, "Gnome", 1, "Exchange", 1, 1, 1, "n", "n^2", "n^2", "c");
    sorts[66] = Sort(66, "Library", 1, "Insertion", 1, 0, 1, "n log n", "n log n", "n^2", "n");
    sorts[67] = Sort(67, "Strand", 1, "Selection", 1, 0, 1, "n", "n^1.5", "n^2", "n");
    sorts[68] = Sort(68, "Patience", 1, "Insertion, Selection", 0, 0, 1, "n", "n log n", "n log n", "n");
    sorts[69] = Sort(69, "Bitonic", 1, "Merge", 0, 1, 1, "n log n log n", "n log n log n", "n log n log n", "c");
    sorts[70] = Sort(70, "Bitonic Network", 1, "Merge", 0, 1, 1, "n log n log n", "n log n log n", "n log n log n", "c");
    sorts[71] = Sort(71, "Pairwise Network", 1, "Merge", 0, 1, 1, "n log n log n", "n log n log n", "n log n log n", "c");

    // OBSCURE NON-COMPARISON SORTS
    sorts[72] = Sort(72, "Spread, 1 bit base", 0, "Buckets", 0, 0, 1, "n", "n log n", "n log n", "n");
    sorts[73] = Sort(73, "Spread, 2 bit base", 0, "Buckets", 0, 0, 1, "n", "n log n", "n log n", "n");
    sorts[74] = Sort(74, "Spread, 4 bit base", 0, "Buckets", 0, 0, 1, "n", "n log n", "n log n", "n");
    sorts[75] = Sort(75, "Spread, 8 bit base", 0, "Buckets", 0, 0, 1, "n", "n log n", "n log n", "n");
    sorts[76] = Sort(76, "Flash, 0.3n buckets", 0, "Buckets", 0, 0, 1, "n", "n", "n^2", "n");
    sorts[77] = Sort(77, "Flash, 0.4n buckets", 0, "Buckets", 0, 0, 1, "n", "n", "n^2", "n");
    sorts[78] = Sort(78, "Flash, 0.5n buckets", 0, "Buckets", 0, 0, 1, "n", "n", "n^2", "n");
    sorts[79] = Sort(79, "Bead", 0, "Buckets", 0, 0, 0, "n", "S", "S", "n^2");
    sorts[80] = Sort(80, "Spaghetti", 0, "Buckets", 1, 0, 0, "n", "n", "n", "n^2");
    sorts[81] = Sort(81, "Ska Sort", 0, "Buckets", 0, 0, 1, "n", "n", "n", "n");
}

void printInfoAllSorts() {
    for (int i = 0; i < numSorts; i++) {
        sorts[i].printInfo();
    }
}

void printInfo(int sort_index) {
    sorts[sort_index].printInfo();
}

void printInfo(string sort_name) {
    bool found = 0;
    for (int i = 0; i < numSorts; i++) {
        if (sorts[i].name == sort_name) {
            sorts[i].printInfo();
            found = 1;
        }
    }
    if (!found) {
        printf("No sorts found with the given name.\n\n");
    }
}

int main(void) {

    srand(0);

    initAllSorts();
    
    n = 100;

    //printInfoAllSorts();

    runAllSorts(1);
}
