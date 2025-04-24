#include <vector>
#include <string>
#include <iostream>
#include <random>

using namespace std;

#include "sort.hpp"

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

    case 15: // pigeonhole
        pigeonhole(v, n);
        break;
    case 16: // counting
        counting(v, n);
        break;
    case 17: // bucket, 4 buckets
        bucket(v, n, 4);
        break;
    case 18: // bucket, 16 buckets
        bucket(v, n, 16);
        break;
    case 19: // bucket, 64 buckets
        bucket(v, n, 64);
        break;
    case 20: // bucket, N buckets
        bucket(v, n, n);
        break;
    case 21: // lsd radix, base 2
        radix_lsd(v, n, 2);
        break;
    case 22: // lsd radix, base 4
        radix_lsd(v, n, 4);
        break;
    case 23: // lsd radix, base 8
        radix_lsd(v, n, 8);
        break;
    case 24: // lsd radix, base 16
        radix_lsd(v, n, 16);
        break;
    case 25: // msd radix, base 2
        radix_msd(v, n, 2);
        break;
    case 26: // msd radix, base 4
        radix_msd(v, n, 4);
        break;
    case 27: // msd radix, base 8
        radix_msd(v, n, 8);
        break;
    case 28: // msd radix, base 16
        radix_msd(v, n, 16);
        break;
    case 29: // in place lsd radix, base 2
        in_place_lsd(v, n, 2);
        break;
    case 30: // in place lsd radix, base 4
        in_place_lsd(v, n, 4);
        break;
    case 31: // in place lsd radix, base 8
        in_place_lsd(v, n, 8);
        break;
    case 32: // in place lsd radix, base 16
        in_place_lsd(v, n, 16);
        break;
    case 33: // in place msd radix, base 2
        in_place_msd(v, n, 2);
        break;
    case 34: // in place msd radix, base 4
        in_place_msd(v, n, 4);
        break;
    case 35: // in place msd radix, base 8
        in_place_msd(v, n, 8);
        break;
    case 36: // in place msd radix, base 16
        in_place_msd(v, n, 16);
        break;

    case 37: // intro
        intro(v, n);
        break;
    case 38: // cycle
        cycle(v, n);
        break;
    case 39: // tim, group size 4
        tim(v, n, 4);
        break;
    case 40: // tim, group size 16
        tim(v, n, 16);
        break;
    case 41: // tim, group size 64
        tim(v, n, 64);
        break;
    case 42: // iterative merge
        iterative_merge(v, n);
        break;
    case 43: // naive in-place merge
        in_place_merge(v, n);
        break;
    case 44: // weave
        weave(v, n);
        break;
    case 45: // iterative weave
        iterative_weave(v, n);
        break;
    case 46: // rotate merge
        rotate_merge(v, n);
        break;
    case 47: // block
        block(v, n);
        break;
    case 48: // iterative block
        iterative_block(v, n);
        break;
    case 49: // wiki
        wiki(v, n);
        break;
    case 50: // grail, block size 4
        grail(v, n, 4);
        break;
    case 51: // grail, block size 16
        grail(v, n, 16);
        break;
    case 52: // grail, block size 64
        grail(v, n, 64);
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
    case 61: // merge-insertion
        merge_insert(v, n);
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
    case 70: // bitonic
        bitonic_network(v, n);
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
    case 79: // spaghetti
        spaghetti(v, n);
        break;
    case 80: // ska
        ska(v, n);
        break;

    case 81: // bead (gravity)
        bead(v, n);
        break;
    case 82: // stooge
        stooge(v, n);
        break;
    case 83: // bogo
        bogo(v, n);
        break;
    }
}

void Sorter::bubble(vector<int>& v, int n) {
    bool swapped;

    for (int i = 0; i < n - 1; i++) {
        swapped = 0;
        for (int j = 0; j < n - i - 1; j++) {
            if (v[j] > v[j + 1]) {
                swap(v[j], v[j + 1]);
                swapped = 1;
            }
        }

        if (!swapped) break;
    }
}

void Sorter::exchange(vector<int>& v, int n)
{
    int i, j;
    for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
            if (v[i] > v[j]) {
                swap(v[i], v[j]);
            }
        }
    }
}

void Sorter::selection(vector<int>& v, int n) {

    for (int i = 0; i < n - 1; ++i) {
        int minI = i;

        for (int j = i + 1; j < n; ++j) {
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
    for (int i = 1; i < n; ++i) {
        int val = v[i];
        int j = i - 1;

        while (j && v[j] > val) {
            v[j + 1] = v[j];
            j--;
        }
        v[j + 1] = val;
    }
}

int Sorter::binary_search(vector<int>& v, int c, int l, int r)
{
    if (r <= l) return (c > v[l]) ? (l + 1) : l;

    int m = (l + r) / 2;

    if (c == v[m]) return m + 1;

    if (c > v[m]) return binary_search(v, c, m + 1, r);
    return binary_search(v, c, l, m - 1);
}

void Sorter::binary_insertion(vector<int>& v, int n)
{
    int i, loc, j, c;

    for (i = 1; i < n; ++i)
    {
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

    vector<int> L(n1), R(n2);

    // Copy data to temp vectors L and R
    for (int i = 0; i < n1; i++) L[i] = v[l + i];
    for (int i = 0; i < n2; i++) R[i] = v[m + 1 + i];

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

void Sorter::merge(vector<int>&v, int n){
    merg(v, 0, n - 1);
}

void Sorter::mer_3(vector<int>&v, int l, int m1, int m2, int r) {
    int s1 = m1 - l + 1;
    int s2 = m2 - m1;
    int s3 = r - m2;

    vector<int> a1(s1), a2(s2), a3(s3);

    // Copy data to temporary arrays
    for (int i = 0; i < s1; i++) {
        a1[i] = v[l + i];
    }
    for (int i = 0; i < s2; i++) {
        a2[i] = v[m1 + 1 + i];
    }
    for (int i = 0; i < s3; i++) {
        a3[i] = v[m2 + 1 + i];
    }

    // Merge sorted subarrays
    int i = 0, j = 0, k = 0, spot = l;
    while (i < s1 || j < s2 || k < s3) {
        int min = INT_MAX;

        // Find smallest
        if (i < s1 && a1[i] < min) {
            min = a1[i++];
        }
        if (j < s2 && a2[j] < min) {
            min = a2[j++];
        }
        if (k < s3 && a3[k] < min) {
            min = a3[k++];
        }

        // Place smallest into original array
        v[spot++] = min;
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
    int n1 = m1 - l;
    int n2 = m2 - m1;
    int n3 = m3 - m2;
    int n4 = r - m3;

    vector<int> L(n1), M(n2), N(n3), R(n4);

    // Copy the subarrays into temporary vectors
    for (int i = 0; i < n1; ++i) L[i] = v[l + i];
    for (int i = 0; i < n2; ++i) M[i] = v[m1 + i];
    for (int i = 0; i < n3; ++i) N[i] = v[m2 + i];
    for (int i = 0; i < n4; ++i) R[i] = v[m3 + i];

    // Merge sorted subarrays
    int h = 0, i = 0, j = 0, k = 0, spot = l;
    while (h < n1 || i < n2 || j < n3 || k < n4) {

        int min = INT_MAX;

        // Find smallest
        if (h < n1 && L[h] < min) {
            min = L[h++];
        }
        if (i < n2 && M[i] < min) {
            min = M[i++];
        }
        if (j < n3 && N[j] < min) {
            min = N[j++];
        }
        if (k < n4 && R[k] < min) {
            min = R[k++];
        }

        // Place smallest into original array
        v[spot++] = min;
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
    merg_4(v, 0, n - 1);
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
    qui_lomuto(v, 0, n - 1);
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
    qui_hoare(v, 0, n - 1);
}

void Sorter::qui_naive(vector<int>& v, int l, int r) {
    if (l >= r) return;

    int pivot = v[r];
    vector<int> L, R;

    for (int i = l; i < r; ++i) {
        if (v[i] <= pivot)
            L.push_back(v[i]);
        else
            R.push_back(v[i]);
    }

    // Copy back sorted values to original vector
    int index = l;
    for (int x:L) v[index++] = x;
    v[index++] = pivot;
    for (int x:R) v[index++] = x;

    int m = l + L.size();
    qui_naive(v, l, m - 1);
    qui_naive(v, m + 1, r);
}

void Sorter::quick_naive(vector<int>& v, int n) {
    qui_naive(v, 0, n - 1);
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
    int gap = n;

    // Initialize swapped as true to make sure that loop runs
    bool swapped = 1;

    while (gap > 1 || swapped) // IDK IF SWAPPED SHOULD BE CHANGED TO !SWAPPED HERE
    {
        // Find next gap
        gap = (gap * 10) / 13;
        // if(gap == 9 || gap == 10)    RULE OF 11, MAYBE ENABLE
        // gap = 11;
        if (gap < 1) gap = 1;

        // Initialize swapped as false so that we can check if swap happened or not
        swapped = 0;

        // Compare all elements with current gap
        for (int i = 0; i < n - gap; i++)
        {
            if (v[i] > v[i + gap])
            {
                swap(v[i], v[i + gap]);
                swapped = 1;
            }
        }
    }
}

void Sorter::shell(vector<int>& v, int n) {
    // Start with a large gap, then reduce it
    for (int gap = n / 2; gap > 0; gap /= 2) {
        // Perform a gapped insertion sort for this gap size
        for (int i = gap; i < n; ++i) {
            int temp = v[i];
            int j;
            for (j = i; j >= gap && v[j - gap] > temp; j -= gap) {
                v[j] = v[j - gap];
            }
            v[j] = temp;
        }
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
    for (int i = 1; i < size(b); ++i) {
        int val = b[i];
        int j = i - 1;
        while (j && b[j] > val) {
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
    float range = (float)(max - min) * 1.001;

    // Put array elements in buckets using their values
    for (int i = 0; i < n; i++) {
        int bi = (float)N * (float)(v[i] - min) / range;
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

void Sorter::radix_count(vector<int>& v, int n, int power, int b){
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
        out[c[(v[i] / power) % b] - 1] = v[i];
        c[(v[i] / power) % b]--;
    }

    v = out;
}

void Sorter::radix_lsd(vector<int>& v, int n, int b){
    // Find the maximum number to know number of digits
    int m = v[0];
    for (int i = 1; i < n; i++)
        if (v[i] > m) m = v[i];

    // Do counting sort for every digit's power
    for (int power = 1; m / power > 0; power *= b) radix_count(v, n, power, b);
}

int Sorter::numDigits(int n, int b){
    int count = 0;
    while (n > 0) n /= b;
    count++;
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

void Sorter::radix_msd(vector<int>& v, int n, int b) {
    // Find the maximum number to know number of digits
    int m = numDigits(v[0], b);
    for (int i = 1; i < n; i++)
        if (numDigits(v[i], b) > m) m = numDigits(v[i], b);

    // Do counting sort for every digit's power
    for (int power = pow(b, m); power > 0; power /= 10) radix_count(v, n, power, b);
}

void Sorter::in_place_lsd(vector<int>& v, int n, int b) {
    if (v.empty()) return;

    // Find max number to determine number of digits
    int maxVal = *max_element(v.begin(), v.end());

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

int Sorter::in_place_msd_digit(int number, int pos, int maxDigits, int b) {
    int divisor = pow(b, maxDigits - pos - 1);
    return (number / divisor) % b;
}

void Sorter::in_place_msd_bucket(vector<int>& v, int start, int end, int pos, int maxDigits, int b) {
    if (end - start <= 1 || pos >= maxDigits) return;

    vector<int> count(b, 0);
    vector<int> offset(b, 0);

    // Count digits at the current position
    for (int i = start; i < end; i++) {
        int d = in_place_msd_digit(v[i], pos, maxDigits, b);
        count[d]++;
    }

    // Get cumulative counts
    offset[0] = start;
    for (int i = 1; i < b; i++) {
        offset[i] = offset[i - 1] + count[i - 1];
    }

    // Rearrange elements in-place
    for (int i = 0; i < b; i++) {
        while (count[i] > 0) {
            int from = offset[i];
            int digit = in_place_msd_digit(v[from], pos, maxDigits, b);
            if (digit == i) {
                offset[i]++;
                count[i]--;
            }
            else {
                int to = offset[digit];
                swap(v[from], v[to]);
                count[digit]--;
                offset[digit]++;
            }
        }
    }

    // Recursively sort each bucket
    for (int i = 0; i < b; i++) {
        int bucketStart = (i == 0) ? start : offset[i - 1];
        int bucketEnd = offset[i];
        in_place_msd_bucket(v, bucketStart, bucketEnd, pos + 1, maxDigits, b);
    }
}

void Sorter::in_place_msd(vector<int>& v, int n, int b) {
    if (n <= 1) return;

    int maxVal = *max_element(v.begin(), v.end());
    int maxDigits = (maxVal == 0) ? 1 : static_cast<int>(log10(maxVal)) + 1;

    in_place_msd_bucket(v, 0, n, 0, maxDigits, b);
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
    // cout << writes << endl ;
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
    // parts left and right array
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
    reverse(v.begin() + l, v.begin() + m);
    reverse(v.begin() + m, v.begin() + r);
    reverse(v.begin() + l, v.begin() + r);
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

    int m = l + (r - l) / 2;
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

bool Sorter::wiki_pass(vector<int>& v, int n) {
    bool swapped = 0;
    for (int i = 0; i < n - 1; ++i) {
        // Swap elements if out of order
        if (v[i] > v[i + 1]) {
            swap(v[i], v[i + 1]);
            swapped = 1;
        }
    }
    return swapped;
}

void Sorter::wiki(vector<int>& v, int n) {
    bool swapped = wiki_pass(v, n);
    while (swapped) {
        swapped = wiki_pass(v, n);
    }
}

void Sorter::grail_insert(vector<int>& v, int l, int r) {
    for (int i = l + 1; i < r; ++i) {
        int key = v[i];
        int j = i - 1;
        while (j >= l && v[j] > key) {
            v[j + 1] = v[j];
            --j;
        }
        v[j + 1] = key;
    }
}

void Sorter::grail_merge(vector<int>& v, int l, int m, int r) {
    int i = l;
    int j = m;
    while (i < j && j < r) {
        if (v[i] <= v[j]) {
            ++i;
        }
        else {
            int temp = v[j];
            for (int k = j; k > i; --k) {
                v[k] = v[k - 1];
            }
            v[i] = temp;
            ++i;
            ++j;
        }
    }
}

void Sorter::grail(vector<int>& v, int n, int k) {

    for (int i = 0; i < n; i += k) {
        grail_insert(v, i, min(i + k, n));
    }

    for (int size = k; size < n; size *= 2) {
        for (int l = 0; l < n; l += 2 * size) {
            int m = min(l + size, n);
            int r = min(l + 2 * size, n);
            grail_merge(v, l, m, r);
        }
    }
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

void Sorter::smooth_build(vector<int>& v, int n, int i) {
    int max = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;

    // Compare with left child
    if (l < n && v[l] > v[max]) {
        max = l;
    }

    // Compare with right child
    if (r < n && v[r] > v[max]) {
        max = r;
    }

    // If the largest is not the root, swap and continue heapifying
    if (max != i) {
        swap(v[i], v[max]);
        smooth_build(v, n, max);
    }
}

void Sorter::smooth(vector<int>& v, int n) {

    make_heap(v.begin(), v.end());  // Builds a heap in O(n) time

    // Extract elements one by one from the heap, similar to Heap Sort
    for (int i = n - 1; i >= 0; --i) {
        swap(v[0], v[i]);
        // "heapify" the remaining unsorted part
        smooth_build(v, i, 0);
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

        // Recursively sort the left and right subarrays
        binary_qui(v, l, pi - 1);  // Left subarray
        binary_qui(v, pi + 1, r); // Right subarray
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

void Sorter::pancake(vector<int>& v, int n) {

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

void Sorter::odd_even(vector<int>& v, int n)
{
    bool isSorted = 0; // Initially array is unsorted

    while (!isSorted)
    {
        isSorted = 1;

        // Perform Bubble sort on odd indexed element
        for (int i = 1; i <= n - 2; i += 2)
        {
            if (v[i] > v[i + 1])
            {
                swap(v[i], v[i + 1]);
                isSorted = 0;
            }
        }

        // Perform Bubble sort on even indexed element
        for (int i = 0; i <= n - 2; i += 2)
        {
            if (v[i] > v[i + 1])
            {
                swap(v[i], v[i + 1]);
                isSorted = 0;
            }
        }
    }

    return;
}

bool Sorter::circle_pass(vector<int>& v, int l, int r) {
    bool swapped = 0;

    // Base case
    if (l >= r) return 0;

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

    int m = (r - l) / 2;
    bool left = circle_pass(v, l, l + m);
    bool right = circle_pass(v, l + m + 1, r);
    return left || right || swapped;
}

void Sorter::circle(vector<int>& v, int n) {
    int pass = 0;

    while (circle_pass(v, pass, n - 1 - pass)) {
        pass++;
    }
}

void Sorter::merge_insert_insert(vector<int>& v, int x) {
    int l = 0, r = size(v);
    while (l < r) {
        int m = l + (r - l) / 2;
        if (x < v[m]) r = m;
        else l = m + 1;
    }
    v.insert(v.begin() + l, x);
}

void Sorter::merge_insert(vector<int>& v, int n) {
    if (n <= 1) return;

    vector<int> out;
    vector<int> temp;

    for (int i = 0; i + 1 < n; i += 2) {
        if (v[i] < v[i + 1]) {
            out.push_back(v[i + 1]);
            temp.push_back(v[i]);
        }
        else {
            out.push_back(v[i]);
            temp.push_back(v[i + 1]);
        }
    }

    if (n % 2 != 0) {
        out.push_back(v[n - 1]);
    }

    merge_insert(out, size(out));

    for (int val : temp) {
        merge_insert_insert(out, val);
    }

    // Copy result back
    v = out;
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

void Sorter::tournament_build(vector<int>& tree, vector<int>&v, int n) {
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

void Sorter::gnome(vector<int>& v, int n)
{
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

    const double epsilon = 1.0;  // Controls extra space (1.0 = 100% extra)
    int s = static_cast<int>((1.0 + epsilon) * n) + 1; // size of vector with gaps

    // Using INT_MAX to represent empty slots (gaps)
    vector<int> temp(s, INT_MAX);

    auto insert = [&](int x) {
        // Binary search for correct insertion index (ignoring gaps)
        int l = 0, r = s;
        while (l < r) {
            int m = (l + r) / 2;
            if (temp[m] == INT_MAX || x < temp[m])
                r = m;
            else
                l = m + 1;
        }

        // Shift right to find actual empty slot
        int pos = r;
        while (pos < s && temp[pos] != INT_MAX)
            ++pos;

        if (pos >= s) {
            return; // out of space while inserting
        }

        // Shift elements to make space
        for (int i = pos; i > r; --i)
            temp[i] = temp[i - 1];
        temp[r] = x;
        };

    // Insert elements one by one
    insert(v[0]);
    for (int i = 1; i < n; ++i) {
        insert(v[i]);
    }

    // Copy back sorted elements (ignoring gaps)
    v.clear();
    for (int x : temp) {
        if (x != INT_MAX) {
            v.push_back(x);
        }
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

vector<int> Sorter::patience_merge(vector<vector<int>>& v)
{
    vector<int> out;

    // Always move smallest (top) to out
    while (1) {

        // Stores the smallest element of the top of the piles
        int min = INT_MAX;

        // Stores index of the smallest element of the top of the piles
        int index = -1;

        // Calculate the smallest element
        for (int i = 0; i < size(v); i++) {

            // If min is greater than top of the current stack
            if (min > v[i][size(v[i]) - 1]) {

                // Update min
                min = v[i][size(v[i]) - 1];

                // Update index
                index = i;
            }
        }

        // Insert the smallest element of the top of the stack
        out.push_back(min);

        // Remove the top element from the current pile
        v[index].pop_back();

        if (v[index].empty()) {

            // Remove current pile from all piles
            v.erase(v.begin() + index);
        }

        // If all the piles are empty
        if (v.size() == 0) break;
    }
    return out;
}

void Sorter::patience(vector<int>& v, int n)
{
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
    v = patience_merge(p);
}

void Sorter::bitonic_merge(vector<int>& v, int low, int n, bool dir){
    if (n > 1)
    {
        int k = n / 2;
        for (int i = low; i < low + k; i++) {
            if (dir == (v[i] > v[i+k]))
                swap(v[i], v[i+k]);
        }

        bitonic_merge(v, low, k, dir);
        bitonic_merge(v, low + k, k, dir);
    }
}

void Sorter::bitonic_sub(vector<int>& v, int low, int n, bool dir){
    if (n > 1)
    {
        int k = n / 2;

        // sort in ascending order
        bitonic_sub(v, low, k, 1);

        // sort in descending order
        bitonic_sub(v, low + k, k, 0);

        // Will merge whole sequence in ascending order since dir = 1
        bitonic_merge(v, low, n, dir);
    }
}

void Sorter::bitonic(vector<int>& v, int n){
    bitonic_sub(v, 0, n, 1);
}

// Recursively merges bitonic sequences
void Sorter::bitonic_network_merge(vector<int>& v, int low, int count, bool dir) {
    if (count <= 1) return;

    int k = count / 2;
    for (int i = low; i < low + k; ++i) {
        if ((dir && v[i] > v[i + k]) || (!dir && v[i] < v[i + k])) {
            swap(v[i], v[i + k]);
        }
    }

    bitonic_network_merge(v, low, k, dir);
    bitonic_network_merge(v, low + k, k, dir);
}

// Main bitonic sort function
void Sorter::bitonic_network_sort(vector<int>& v, int low, int count, bool dir) {
    if (count <= 1) return;

    int k = count / 2;

    // Sort first half ascending and second half descending
    bitonic_network_sort(v, low, k, 1);
    bitonic_network_sort(v, low + k, k, 0);

    // Merge whole sequence in desired direction
    bitonic_network_merge(v, low, count, dir);
}

void Sorter::bitonic_network(vector<int>& v, int n) {

    int power = 1;
    while (power < n) power *= 2;
    for (int i = n; i < power; ++i) v.push_back(INT_MAX); // Padding

    bitonic_network_sort(v, 0, size(v), 1);

    // Remove padding
    v.resize(n);
}

void Sorter::pairwise_network_swap(vector<int>& v, int i, int j) {
    if (v[i] > v[j]) {
        swap(v[i], v[j]);
    }
}

void Sorter::pairwise_network_merge(vector<int>& v, int l, int r) {
    int m = (l + r) / 2;
    if (r - l <= 1) return;

    for (int i = l; i < m; ++i) {
        pairwise_network_swap(v, i, i + (m - l));
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
    // Pad size to power of 2
    int power = 1;
    while (power < n) power *= 2;
    for (int i = n; i < power; ++i) v.push_back(INT_MAX); // Padding with large values

    pairwise_network_sort(v, 0, v.size());

    // Remove padding
    v.resize(n);
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
    for (int i = 0; i < bucketCount; ++i) {
        if (!buckets[i].empty()) {
            if (size(buckets[i]) > bucketCount) {
                // Recursively sort large buckets
                for (int val: buckets[i]) v[index++] = val;
                spread_sort(v, index - size(buckets[i]), index - 1, shift - b, b);
            }
            else {
                // Sort for small buckets
                sort(buckets[i].begin(), buckets[i].end());
                for (int val: buckets[i]) v[index++] = val;
            }
        }
    }
}

void Sorter::spread(vector<int>& v, int n, int b) {
    if (n == 0) return;

    int max = *max_element(v.begin(), v.end());
    int shift = spread_digits(max);

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

void Sorter::spaghetti(vector<int>& v, int n) {

    int max = *max_element(v.begin(), v.end());

    vector<int> rods(max + 1, 0);

    for (int i = 0; i < n; ++i) ++rods[v[i]];

    // Collect spaghetti in descending count order
    int j = 0;
    for (int i = 0; i <= max; ++i) {
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

void Sorter::bead(vector<int>& v, int n) {
    if (n == 0) return;

    int maxVal = *max_element(v.begin(), v.end());
    if (maxVal == 0) return;

    // Create a grid of "beads"
    vector<vector<bool>> beads(n, vector<bool>(maxVal, false));

    // Place beads on the rods
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < v[i]; ++j) {
            beads[i][j] = true;
        }
    }

    // Simulate gravity
    for (int j = 0; j < maxVal; ++j) {
        int count = 0;
        // Count beads in each column
        for (int i = 0; i < n; ++i)
            if (beads[i][j])
                ++count;

        // Drop beads to the bottom
        for (int i = n - 1; i >= 0; --i) {
            beads[i][j] = (count-- > 0);
        }
    }

    // Read off the sorted values
    for (int i = 0; i < n; ++i) {
        int count = 0;
        for (int j = 0; j < maxVal; ++j)
            if (beads[i][j])
                ++count;
        v[i] = count;
    }
}

void Sorter::stooge_sub(vector<int>& v, int l, int r)
{

    if (v[l] > v[r]) swap(v[r], v[l]);

    if (r - l + 1 > 2)
    {
        int m = (r - l + 1) / 3;
        stooge_sub(v, l, r - m);
        stooge_sub(v, l + m, r);
        stooge_sub(v, l, r - m);
    }
}

void Sorter::stooge(vector<int>& v, int n)
{
    stooge_sub(v, 0, n - 1);
}

void Sorter::bogo(vector<int>& v, int n) {
    int N; bool sorted;
    while (1) {
        N = n;
        sorted = 1;

        // if array is not sorted then shuffle the array again
        while (--n > 0)
            if (v[n] < v[n - 1]) {
                sorted = 0;
                break;
            }
        if (sorted) return;

        for (int i = 0; i < n; i++) swap(v[i], v[rand() % n]);
    }
}