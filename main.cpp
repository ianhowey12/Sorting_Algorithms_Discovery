#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <chrono>

using namespace std;

#include "sort.hpp"
#include "descriptions.hpp"

#define numSorts 84
#define numAlgorithms 59

vector<Sort> sorts(numSorts);

// vector to be sorted and its size
vector<int> v;
int n;

void fillVector() {
    while (size(v) > 0) v.pop_back();

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

double runSort(bool print, int i) {
    Sorter sorter = Sorter();
    fillVector();

    auto start = chrono::high_resolution_clock::now();

    sorter.sortVector(i, v, n);

    auto end = chrono::high_resolution_clock::now();

    string status = ") failed. ";
    bool success = 0;

    if (checkForSorted()) {
        status = ") succeeded. ";
        success = 1;
    }

    chrono::nanoseconds diff = end - start;
    double timer = (double)diff.count() / 1000000000.0;

    if (print) {
        cout << "Sort #" << i + 1 << " (" << sorts[i].name << status << "Time taken: " << timer << " seconds.\n\n";
    }
    if (!success) timer *= -1.0;
    return timer;
}

// run every sorting algorithm except the slow ones
void runAllSorts(bool print) {

    if (n < 1) {
        cout << "Number of vector items (n) must be at least 1.\n\n";
        return;
    }
    if (n > 10000000) {
        cout << "Number of vector items (n) must be at most 10000000.\n\n";
        return;
    }

    Sorter sorter = Sorter();

    for (int i = 0; i < numSorts - 3; i++) {
        runSort(print, i);
    }
}

// run the slow sorting algorithms
void runSlowSorts(bool print) {

    if (n < 1) {
        cout << "Number of vector items (n) must be at least 1.\n\n";
        return;
    }
    if (n > 100) {
        cout << "Number of vector items (n) must be at most 100 for slow sorts.\n\n";
        return;
    }

    for (int i = numSorts - 3; i < numSorts; i++) {
        runSort(print, i);
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

    // POPULAR NON-COMPARISON SORTS
    sorts[15] = Sort(15, "Pigeonhole", 0, "Buckets", 1, 0, 0, "", "", "", "");
    sorts[16] = Sort(16, "Counting", 0, "Buckets", 1, 0, 0, "n+b^d", "n+b^d", "n+b^d", "n+b^d");
    sorts[17] = Sort(17, "Bucket, 4 Buckets", 0, "Buckets", 1, 0, 1, "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b");
    sorts[18] = Sort(18, "Bucket, 16 Buckets", 0, "Buckets", 1, 0, 1, "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b");
    sorts[19] = Sort(19, "Bucket, 64 Buckets", 0, "Buckets", 1, 0, 1, "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b+(n^2)/b", "n+b");
    sorts[20] = Sort(20, "Bucket, N Buckets", 0, "Buckets", 1, 0, 1, "3n", "3n", "3n", "2n");
    sorts[21] = Sort(22, "LSD Radix, Base 2", 0, "Buckets", 1, 0, 1, "nd+bd", "nd+bd", "nd+bd", "n+b");
    sorts[22] = Sort(23, "LSD Radix, Base 4", 0, "Buckets", 1, 0, 1, "nd+bd", "nd+bd", "nd+bd", "n+b");
    sorts[23] = Sort(24, "LSD Radix, Base 8", 0, "Buckets", 1, 0, 1, "nd+bd", "nd+bd", "nd+bd", "n+b");
    sorts[24] = Sort(25, "LSD Radix, Base 16", 0, "Buckets", 1, 0, 1, "nd+bd", "nd+bd", "nd+bd", "n+b");
    sorts[25] = Sort(26, "MSD Radix, Base 2", 0, "Buckets", 1, 0, 1, "ndb", "ndb", "ndb", "n+b+d");
    sorts[26] = Sort(27, "MSD Radix, Base 4", 0, "Buckets", 1, 0, 1, "ndb", "ndb", "ndb", "n+b+d");
    sorts[27] = Sort(28, "MSD Radix, Base 8", 0, "Buckets", 1, 0, 1, "ndb", "ndb", "ndb", "n+b+d");
    sorts[28] = Sort(29, "MSD Radix, Base 16", 0, "Buckets", 1, 0, 1, "ndb", "ndb", "ndb", "n+b+d");
    sorts[29] = Sort(30, "In-Place LSD Radix, Base 2", 0, "Buckets", 1, 1, 1, "dn^2", "dn^2", "dn^2", "b");
    sorts[30] = Sort(31, "In-Place LSD Radix, Base 4", 0, "Buckets", 1, 1, 1, "dn^2", "dn^2", "dn^2", "b");
    sorts[31] = Sort(32, "In-Place LSD Radix, Base 8", 0, "Buckets", 1, 1, 1, "dn^2", "dn^2", "dn^2", "b");
    sorts[32] = Sort(33, "In-Place LSD Radix, Base 16", 0, "Buckets", 1, 1, 1, "dn^2", "dn^2", "dn^2", "b");
    sorts[33] = Sort(34, "In-Place MSD Radix, Base 2", 0, "Buckets", 1, 1, 1, "ndb", "ndb", "ndb", "b+d");
    sorts[34] = Sort(35, "In-Place MSD Radix, Base 4", 0, "Buckets", 1, 1, 1, "ndb", "ndb", "ndb", "b+d");
    sorts[35] = Sort(36, "In-Place MSD Radix, Base 8", 0, "Buckets", 1, 1, 1, "ndb", "ndb", "ndb", "b+d");
    sorts[36] = Sort(37, "In-Place MSD Radix, Base 16", 0, "Buckets", 1, 1, 1, "ndb", "ndb", "ndb", "b+d");

    // OBSCURE COMPARISON SORTS
    sorts[37] = Sort(38, "Intro", 1, "Partitioning, Selection", 0, 0, 1, "n log n", "n log n", "n log n", "log n");
    sorts[38] = Sort(39, "Cycle", 1, "Selection", 0, 1, 1, "n^2", "n^2", "n^2", "c");
    sorts[39] = Sort(40, "Tim, 4", 1, "Merge", 1, 0, 1, "n", "n log n", "n log n", "n");
    sorts[40] = Sort(41, "Tim, 16", 1, "Merge", 1, 0, 1, "n", "n log n", "n log n", "n");
    sorts[41] = Sort(42, "Tim, 64", 1, "Merge", 1, 0, 1, "n", "n log n", "n log n", "n");
    sorts[42] = Sort(43, "Iterative Merge", 1, "Merge", 1, 0, 1, "n log n", "n log n", "n log n", "n");
    sorts[43] = Sort(44, "Naive In-Place Merge", 1, "Merge", 1, 1, 1, "n^2", "n^2", "n^2", "c");
    sorts[44] = Sort(45, "Weave", 1, "Merge", 0, 1, 1, "n^2", "n^2", "n^2", "c");
    sorts[45] = Sort(45, "Iterative Weave", 1, "Merge", 0, 1, 1, "n^2", "n^2", "n^2", "c");
    sorts[46] = Sort(46, "Rotate Merge", 1, "Merge", 1, 1, 1, "n log n", "n^2", "n^2", "c");
    sorts[47] = Sort(47, "Block", 1, "Merge", 1, 1, 1, "n log n", "n log n", "n log n", "c");
    sorts[48] = Sort(48, "Iterative Block", 1, "Merge", 1, 1, 1, "n log n", "n log n", "n log n", "c");
    sorts[49] = Sort(49, "Wiki", 1, "Exchange", 1, 1, 1, "n", "n^2", "n^2", "c");
    sorts[50] = Sort(50, "Grail, 4", 1, "Merge", 0, 0, 1, "n log n", "n log n", "n^2", "n");
    sorts[51] = Sort(51, "Grail, 16", 1, "Merge", 0, 0, 1, "n log n", "n log n", "n^2", "n");
    sorts[52] = Sort(52, "Grail, 64", 1, "Merge", 0, 0, 1, "n log n", "n log n", "n^2", "n");
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
    sorts[66] = Sort(66, "Library", 1, "Insertion", 0, 0, 1, "n", "n log n", "n^2", "n");
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
    sorts[79] = Sort(79, "Spaghetti", 0, "Buckets", 1, 0, 0, "n", "n", "n", "n^2");
    sorts[80] = Sort(80, "Ska", 0, "Buckets", 0, 0, 1, "n", "n", "n", "n");

    sorts[81] = Sort(81, "Bead", 0, "Buckets", 0, 0, 0, "n^2", "n^2", "n^2", "n^2");
    sorts[82] = Sort(82, "Stooge", 1, "Exchange", 0, 0, 1, "n", "n^2.7", "n^2.7", "n");
    sorts[83] = Sort(83, "Bogo", 1, "Exchange", 0, 1, 1, "n", "n!n", "infinite", "c");
}

void printInfoAllSorts() {
    for (int i = 0; i < numSorts; i++) {
        sorts[i].printInfo();
    }
    cout << "\n";
}

void printInfo(int sort_index) {
    if (sort_index < 0) {
        cout << "Info index must be at least 0.\n\n";
        return;
    }
    if (sort_index >= numSorts) {
        cout << "Info index must be at most " << numSorts << ".\n\n";
        return;
    }
    sorts[sort_index].printInfo();
    cout << "\n";
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
        cout << "No sorts found with the given name.\n\n";
    }
}

void printDescriptionAllSorts() {
    for (int i = 0; i < numAlgorithms; i++) {
        cout << algorithmNames[i] << ": " << algorithmDescriptions[i] << "\n\n";
    }
    cout << "\n\n";
}

void printDescription(int sort_index) {
    if (sort_index < 0) {
        cout << "Description index must be at least 0.\n\n";
        return;
    }
    if (sort_index >= numAlgorithms) {
        cout << "Description index must be at most " << numAlgorithms << ".\n\n";
        return;
    }
    cout << algorithmNames[sort_index] << ": " << algorithmDescriptions[sort_index];
    cout << "\n\n";
}

void printDescription(string sort_name) {
    bool found = 0;
    for (int i = 0; i < numAlgorithms; i++) {
        if (algorithmNames[i] == sort_name) {
            cout << algorithmNames[i] << ": " << algorithmDescriptions[i];
            cout << "\n\n";
        }
    }
    if (!found) {
        cout << "No sorts found with the given name.\n\n";
    }
}

// define your own custom sort here
double runCustomSort(bool print) {
    Sorter sorter = Sorter();
    fillVector();

    auto start = chrono::high_resolution_clock::now();

    // change this line to define your own custom sort
    sorter.radix_lsd(v, n, 4);

    auto end = chrono::high_resolution_clock::now();

    string status = " failed. ";
    bool success = 0;

    if (checkForSorted()) {
        status = " succeeded. ";
        success = 1;
    }

    chrono::nanoseconds diff = end - start;
    double timer = (double)diff.count() / 1000000000.0;

    if (print) {
        cout << "Custom sort" << status << "Time taken: " << timer << " seconds.\n\n";
    }
    if (!success) timer *= -1.0;
    return timer;
}

void getBarchart(int startingN, int endingN, int numBars, int numTrialsPerBar, int chartHeight) {
    if (startingN < 1) {
        cout << "The starting N value is " << startingN << " and it must be at least 1.\n";
        return;
    }
    if (endingN < 1) {
        cout << "The ending N value is " << endingN << " and it must be at least 1.\n";
        return;
    }
    if (startingN > endingN) {
        cout << "The starting N value cannot be larger than the ending N value.\n";
        return;
    }
    if (numBars < 1) {
        cout << "The number of bars in the chart is " << numBars << " and it must be at least 1.\n";
        return;
    }
    if (numBars > 200) {
        cout << "The number of bars in the chart is " << numBars << " and it must be at most 200.\n";
        return;
    }

    vector<double> times(numBars, 0.0);
    vector<int> Ns(numBars, 0);

    auto chart_start = chrono::high_resolution_clock::now();

    for (int i = 0; i < numBars; i++) {
        Ns[i] = 0.01 + ((double)i / (double)numBars) * (double)(endingN - startingN) + (double)startingN;
        n = Ns[i];
        for (int j = 0; j < numTrialsPerBar; j++) {
            times[i] += abs(runCustomSort(0));
        }
        times[i] /= (double)numTrialsPerBar;
    }

    auto chart_end = chrono::high_resolution_clock::now();

    chrono::nanoseconds diff = chart_end - chart_start;
    double timer = (double)diff.count() / 1000000000.0;

    vector<int> scaledTimes(numBars, 0);
    double maxTime = -1.0;
    double minTime = (double)INT_MAX;
    int maxTimeIndex = 0;
    int minTimeIndex = 0;
    for (int i = 0; i < numBars; i++) {
        if (times[i] > maxTime) {
            maxTime = times[i];
            maxTimeIndex = i;
        }
        if (times[i] < minTime) {
            minTime = times[i];
            minTimeIndex = i;
        }
    }
    for (int i = 0; i < numBars; i++) {
        scaledTimes[i] = 0.01 + (double)chartHeight * times[i] / maxTime;
    }
    
    for (int i = 0; i < chartHeight; i++) {
        cout << "|";
        for (int j = 0; j < numBars; j++) {
            if (scaledTimes[j] >= chartHeight - i) {
                cout << "#";
            }
            else {
                cout << " ";
            }
        }
        cout << "|\n";
    }
    cout << "-";
    for (int i = 0; i < numBars; i++) {
        cout << "-";
    }
    cout << "-\n\n";

    cout << "Displaying the time distribution of the custom sort for " << startingN << " <= n <= " << endingN << " over " << numBars << " bars, with " << numTrialsPerBar << " trials per bar.\n";
    cout << "The difference in values of n between consecutive bars is " << (double)(endingN - startingN) / (double)numBars << ".\n";
    cout << "The shortest average time (" << minTime << " seconds) was achieved at n = " << Ns[minTimeIndex] << ".\n";
    cout << "The longest average time (" << maxTime << " seconds) was achieved at n = " << Ns[maxTimeIndex] << ".\n";
    cout << "This entire procedure involved " << numBars * numTrialsPerBar << " sorting trials and took " << timer << " seconds to complete.\n\n\n";
}

int main(void) {

    srand(0);

    initAllSorts();

    n = 5000;

    //printInfoAllSorts();

    //printDescriptionAllSorts();

    auto full_start = chrono::high_resolution_clock::now();

    runSort(1, 12);
    //runAllSorts(1);
    //runCustomSort(1);

    auto full_end = chrono::high_resolution_clock::now();

    chrono::nanoseconds diff = full_end - full_start;
    double timer = (double)diff.count() / 1000000000.0;
    cout << "\nAll sorts have been finished in " << timer << " seconds.\n\n\n";


    getBarchart(100, 10000, 100, 1, 30);
}