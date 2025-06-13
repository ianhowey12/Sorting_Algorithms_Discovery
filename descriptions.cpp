#include <vector>
#include <string>
#include <iostream>

using namespace std;

#include "sort.hpp"

static vector<string> descriptions = {


	// Bubble
	"Move left to right, swapping pairs of elements that are out of order. Repeat.",

	// Exchange
	"For every element from left to right, compare it to every element to the right of it and swap appropriately.",

	// Selection
	"Find the minimum element in the unsorted part of the array. Swap it with the first element. Repeat.",

	// Double Selection
	"Find the minimum element and maximum element in the unsorted part of the array. Swap them with the first element and last element. Repeat.",

	// Insertion
	"Starting from the second element and increasing, move every element to the left right one position until reaching a value less than or equal to the original value. Place the original element here. Repeat.",

	// Binary Insertion
	"Starting from the second element and increasing, use binary search to find the insertion spot. Move every element to the right over. Place the original element here. Repeat.",

	// Merge
	"Starting from smaller subsections and working up to the full array, divide the array into two parts. Concatenate the parts in sorted order.",

	// Three-Way Merge
	"Starting from smaller subsections and working up to the full array, divide the array into three parts. Concatenate the parts in sorted order.",

	// Four-Way Merge
	"Starting from smaller subsections and working up to the full array, divide the array into four parts. Concatenate the parts in sorted order.",

	// Quick, Lomuto Partition
	"Starting from smaller subsections and working up to the full array, partition the array into two parts using the Lomuto partition method. Place the small elements on the left side and large elements on the right side.",

	// Quick, Hoare Partition
	"Starting from smaller subsections and working up to the full array, partition the array into two parts using the Hoare partition method. Place the small elements on the left side and large elements on the right side.",

	// Quick, Naive Partition
	"Starting from smaller subsections and working up to the full array, partition the array into two parts using a naive partition method. Place the small elements on the left side and large elements on the right side.",

	// Heap
	"Build a heap from every element. Find the largest element by swapping elements and rebuilding the heap. Move the largest element to the end of the heap and repeat.",

	// Comb
	"Compare all elements with elements in a gap. Swap appropriately. Reduce the gap size.",

	// Shell
	"Move elements in a gap by a multiple of the gap size as an insertion method. Reduce the gap size.",

	// Bogo
	"If the array is not sorted, swap every element with a random element. Repeat.",


	// Pigeonhole
	"Assign holes for every value from the minimum to the maximum value. Insert every element in its value's hole. Concatenate the holes.",

	// Bucket
	"Place every element in a bucket corresponding to a range of numbers containing its value. Sort the buckets using Insertion. Concatenate the buckets.",

	// Counting
	"Count the number of elements with each value from the minimum to the maximum using separate counts. Compute the cumulative count of every count. Add values to new output.",

	// LSD Radix
	"Starting from the least significant digit, put every element in a bucket corresponding to the value of a digit in an element. Concatenate the buckets.",

	// MSD Radix
	"Starting from the most significant digit, put every element in a bucket corresponding to the value of a digit in an element. Concatenate the buckets.",

	// In-Place LSD Radix
	"Starting from the least significant digit, put every element in a bucket corresponding to the value of a digit in an element. Sort the buckets recursively while maintaining relative order. Concatenate the buckets.",

	// In-Place MSD Radix
	"Starting from the most significant digit, put every element in a bucket corresponding to the value of a digit in an element. Sort the buckets recursively while maintaining relative order. Concatenate the buckets.",


	// Intro
	"Partition the array recursively. Sort small sections using insertion. Build a heap, swap elements in the heap and rebuild the heap for other sections that have been subdivided a set number of times.",

	// Cycle
	"For every element, find the number of smaller elements in the part of the array to the right of the element. Sort the rest of this element's cycle in the array.",

	// Tim
	"Sort the array's subarrays of a predefined size using insertion. Merge the subarrays together until reaching the original array size.",

	// Iterative Merge
	"Starting from smaller subsections and working up to the full array iteratively, divide the array into two parts. Concatenate the parts in sorted order.",

	// Naive In-Place Merge
	"Starting from smaller subsections and working up to the full array using two pointers, divide the array into two parts. Concatenate the parts in sorted order.",

	// Weave
	"Starting from smaller subsections and working up to the full array, divide the array into two parts. Sort the first part ascending and the second part descending.",

	// Iterative Weave
	"Starting from smaller subsections and working up to the full array, divide the array into two parts. Sort the first part ascending and the second part descending iteratively.",

	// Rotate Merge
	"Starting from smaller subsections and working up to the full array, divide the array into two parts. Rotate the right part into order.",

	// Block
	"Starting from smaller subsections and working up to the full array, find the intended position of each element. Rotate each subsection into place.",

	// Wiki
	"Swap every element with the next element if needed. Repeat until no swaps occur on one pass.",

	// Grail
	"Starting from smaller subsections and working up to the full array, move every element to the left or right of the middle element by comparing with the middle element. Combine the left and right arrays.",

	// Stooge
	"Starting from smaller subsections and working up to the full array, swap the left and right elements if needed.",

	// Weak Heap
	"Build a weak heap from every element. Find the largest element by swapping elements and rebuilding the weak heap. Move the largest element to the end of the weak heap and repeat.",

	// Smooth
	"Build a smooth heap from every element. Find the largest element by swapping elements and rebuilding the smooth heap. Move the largest element to the end of the smooth heap and repeat.",

	// Poplar Heap
	"Build a poplar heap from every element. Find the largest element by swapping elements and rebuilding the poplar heap. Move the largest element to the end of the poplar heap and repeat.",

	// Binary Quick
	"Starting from smaller subsections and working up to the full array, move all elements smaller than the last element to the left of it. Move the pivot to its correct location.",

	// Pancake
	"Find the largest unsorted element. Reverse the part of the array up to this element. Reverse the whole unsorted part of the array. Repeat.",

	// Cocktail Shaker
	"Move left to right, swapping pairs of elements that are out of order. Move right to left, swapping pairs of elements that are out of order. Repeat.",

	// Odd-Even
	"Swap all odd-indexed adjacent pairs of elements. Swap all even-indexed adjacent pairs of elements. Repeat until sorted.",

	// Circle
	"Swap the first and last element if they are out of order. Move the first and last indices inward and repeat. Repeat this entire pass starting one spot inward until no swaps are made on one pass.",

	// Merge-Insertion
	"Starting from smaller subsections and working up to the full array, sort small subsections using insertion and merge them.",

	// Tree
	"Build a binary tree containing every element. Find each next smallest element by traversing the tree.",

	// Tournament
	"Build a binary tree containing every element. Find each next smallest element by traversing the tree and remove it.",

	// Stable Tournament
	"Build a binary tree containing every element. Find each next smallest element by traversing the tree and remove it. Use stored indices to order duplicate elements.",

	// Gnome
	"Increment an index until the value at that index and the previous value are out of order. Swap them and decrement the index. Repeat until the index reaches the end of the array.",

	// Library
	"For every element, use binary search to find its intended location. Shift elements to create a gap.",

	// Strand
	"Extract an increasing subsequence from the input array. Use two iterators to merge this strand with the output array. Repeat until the input array is empty.",

	// Patience
	"Insert each element into a pile with a top element greater than or equal to this element. Find the smallest element from the tops of the piles and add it to output until all piles are empty.",

	// Bitonic
	"Divide input in half. Sort both halves increasing and decreasing. Merge.",

	// Odd_Even Merge
	"Divide input in half. Sort both halves by comparing and swapping pairs. Merge using a network.",

	// Pairwise Network
	"Divide input in half. Sort both halves by comparing and swapping pairs. Merge using a network.",


	// Spread
	"Split all elements into buckets based on their first digits. Repeat within each bucket until reaching single-value buckets. For small buckets, use a comparison sort. Concatenate the buckets.",

	// Flash
	"Create buckets corresponding to element value ranges. Insert all elements into buckets. Permute elements to approximately correct positions. Move elements back to original array. Use insertion to sort array.",

	// Bead (Gravity)
	"Create an array of beads with the number of beads in each row corresponding to the value of each element. Move all beads downwards. Add the beads in each row and fill each element of the original array with this sum.",

	// Spaghetti
	"Create a count of elements with each value from the minimum to the maximum value in the array. Remove the values from each count bucket in order and replace the original array values with them.",

	// Ska
	"For every digit within the element values, make buckets. Increment the buckets corresponding to this digit in every value. Compute the cumulative count for each bucket. Fill a temporary array with the bucket values at locations defined by the counts. Move all elements to the original array."

};
