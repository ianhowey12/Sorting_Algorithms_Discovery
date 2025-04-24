#include <vector>
#include <string>
#include <iostream>

using namespace std;

#include "sort.hpp"

//Sort::Sort(){}

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

string boolString(bool b) {
    if (b) return "Yes";
    return "No";
}

void Sort::printInfo() {
    cout << "Sort #" << index + 1 << ": " << name << "\nIs a comparison sort: " << boolString(isComparison) << "\nAlgorithm Type: " << methodType << "\nPreserves order of duplicate elements (stable): "
    << boolString(preservesOrder) << "\nSorts elements in-place (without additional storage): " << boolString(inPlace) << "\nCan be modified to sort decimals: " << boolString(canSortDecimals) << "\nBest-case time complexity: " << bestTimeComplexity
    << "\nAverage-case time complexity: " << averageTimeComplexity << "\nWorst-case time complexity: " << worstTimeComplexity << "\nSpace complexity: " << spaceComplexity << "\n\n";
}