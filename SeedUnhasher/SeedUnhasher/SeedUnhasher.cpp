#include <iostream>
#include <string>
#include <climits>
using namespace std;

int MAX_INT = 0x7FFFFFFF;
long long OVERFLOW_CANCEL = 2LL * (long long)MAX_INT + 2LL;

int minChar = 0;
int maxChar = 0;

int hashCode(string s);
long long pow(int base, int exponent);
void printString(int* string, int size);
void recursiveSearch(int maxDepth, long long maxHash, int currentDepth, long long currentHash, int* currentString);
void startSearch(int hash, int minSize, int maxSize);

int main() {
	while(true) {
		cout << "Enter a seed to unhash : ";
		long long gameSeed;
		cin >> gameSeed;

		if (gameSeed > MAX_INT) {
			cout << "The seed you entered was not generated via a string.\n\n";
			continue;
		}

		//Sets the range of characters to find. In this case, it is [ ] through [~].
		minChar = 32;
		maxChar = 126;

		//Looks for strings with the same hash with size 1 through 8.
		startSearch(gameSeed, 1, 8);
	}

    return 0;
}

int hashCode(string s) {
	int hash = 0;
	int size = s.size();

	for (int i = 0; i < size; ++i) {
		hash = (hash * 31) + s.at(i);
	}

	return hash;
}

long long pow(int base, int exponent) {
	long long value = 1;
	for (int i = 0; i < exponent; ++i) {
		value *= base;
	}
	return value;
}

void printString(int* string, int size) {
	cout << '[';
	for(int i = 0; i < size; ++i) {
		cout << (char)(string[i]);
	}
	cout << ']' << endl;
}

int calculateOverflowFrenquency(int size) {
	long long highestHash = 0;
	int overflowCount = 0;

	for (int i = 0; i < size; ++i) {
		highestHash = (highestHash * 31) + maxChar;
	}

	while(highestHash > MAX_INT) {
		highestHash -= MAX_INT;
		++overflowCount;
	}

	return overflowCount;
}

void startSearch(int hash, int minSize, int maxSize) {
	for(int size = minSize; size <= maxSize; ++size) {
		int overflowFrequency = calculateOverflowFrenquency(size);
		long long searchHash = hash;
		for(long long i = 0; i <= overflowFrequency; ++i) {
			recursiveSearch(size - 1, searchHash + OVERFLOW_CANCEL * i, 0, 0, new int[size - 1]);
		}
	}
}

void recursiveSearch(int maxDepth, long long maxHash, int currentDepth, long long currentHash, int* currentString) {
	if(currentDepth == maxDepth) {
		int lastChar = maxHash - currentHash;
		if(lastChar < minChar || lastChar > maxChar)return;
		currentString[currentDepth] = lastChar;
		printString(currentString, maxDepth + 1);
		return;
	}

	int currentExponent = maxDepth - currentDepth;
	long long currentMultiplier = pow(31, currentExponent);

	long long minAddEnd = 0;
	long long maxAddEnd = 0;

	for(int i = currentExponent - 1; i > -1; --i) {
		long long multiplier = pow(31, i);
		minAddEnd += (long long)minChar * multiplier;
		maxAddEnd += (long long)maxChar * multiplier;
	}

	int currentMinChar = (maxHash - currentHash - maxAddEnd) / currentMultiplier;
	int currentMaxChar = (maxHash - currentHash - minAddEnd) / currentMultiplier;

	if(currentMinChar < minChar || currentMaxChar > maxChar)return;

	for(int i = currentMinChar; i <= currentMaxChar; ++i) {
		currentString[currentDepth] = i;
		recursiveSearch(maxDepth, maxHash, currentDepth + 1, currentHash + ((long long)i * currentMultiplier), currentString);
	}
}

