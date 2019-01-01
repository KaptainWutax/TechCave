#pragma once
#include "Math.h"

//Standard librairies
#include <utility>
#include <vector>
using namespace std;

class PatternSearch {
private:
	int* pattern;
	int patternHeight;
	int patternWidth;

	int *matchTable;

	//These fields are used for debugging, not for the actual search.
	double *chanceTable;
	double averageSkip;

	Math* math;

	void startSearchJava(long long gameSeed, vector<pair<int, int>> boundaries, bool searchForAll, bool showDebugInfo);
	void startSearchBedrock(long long gameSeed, vector<pair<int, int>> boundaries, bool searchForAll, bool showDebugInfo);

public:
	enum Platform {
		JAVA_EDITION,
		BEDROCK_EDITION
	};

	PatternSearch(int h, int w, vector<int> p);
	~PatternSearch();

	void startSearch(Platform platform, long long gameSeed, vector<pair<int, int>> boundaries, bool searchForAll, bool showDebugInfo);

	void preComputeMatchTable();

	//These methods are used for debugging, not for the actual search.
	void preComputeChanceTable();
	void preComputeAverageSkip();

	void printPattern();
	void printMatchTable();
	void printChanceTable();
	void printAverageSkip();
};

