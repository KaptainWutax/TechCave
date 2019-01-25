#include "PatternSearch.h"
#include "Math.h"

//Standard librairies
#include <chrono>
#include <iostream>
#include <thread>
#include <utility>
#include <vector>
using namespace std;

PatternSearch::PatternSearch(int h, int w, vector<int> p) : patternHeight(h), patternWidth(w) {
	if(p.size() != patternHeight * patternWidth || patternHeight * patternWidth == 0)delete this;

	this->pattern = new int[p.size()];
	for(int patternIndex = 0; patternIndex < p.size(); ++patternIndex) {
		this->pattern[patternIndex] = p.at(patternIndex);
	}

	this->math = new Math(this->patternHeight);

	this->preComputeMatchTable();
	this->preComputeChanceTable();
	this->preComputeAverageSkip();
}

PatternSearch::~PatternSearch() {
}

//--------------={SEARCH STARTS HERE}=--------------//
void PatternSearch::startSearch(Platform platform, long long gameSeed, vector<pair<int, int>> boundaries, bool searchForAll, bool showDebugInfo) {
	if(platform == this->JAVA_EDITION) {
		startSearchJava(gameSeed, boundaries, searchForAll, showDebugInfo);
	} else if(platform == this->BEDROCK_EDITION) {
		startSearchBedrock(gameSeed, boundaries, searchForAll, showDebugInfo);
	}
}

void PatternSearch::startSearchJava(long long gameSeed, vector<pair<int, int>> boundaries, bool searchForAll, bool showDebugInfo) {	
	if(showDebugInfo) {
		this->printPattern();
		this->printMatchTable();
		this->printChanceTable();
		this->printAverageSkip();
	}

	int minX = boundaries.at(0).first;
	int maxX = boundaries.at(0).first;
	int minZ = boundaries.at(0).second;
	int maxZ = boundaries.at(0).second;

	int x = boundaries.at(1).first;
	int z = boundaries.at(1).second;
	if(minX > x)minX = x;
	if(maxX < x)maxX = x;
	if(minZ > z)minZ = z;
	if(maxZ < z)maxZ = z;

	int nmOfPatternsFound = 0;
	auto startTime = chrono::system_clock::now();

	int Xskip = 1;
	for(int chunkZ = minZ; chunkZ != maxZ + 1; ++chunkZ) {
		for(int chunkX = minX; chunkX < maxX + 1; chunkX += Xskip) {
			int startSum = 0;
			bool isNotMatch = false;

			for(int i = this->patternHeight - 1; i != -1; --i) {
				//-----={Start check if slime chunk}=-----//
				int x = chunkX + this->patternWidth - 1;
				int z = chunkZ + i;

				long long seed =
					gameSeed +
					(long long)(x * x * 4987142) +
					(long long)(x * 5947611) +
					((long long)(z * z) * (long long)4392871) +
					(long long)(z * 0x5F24F) ^ (long long)987234911;

				seed = (seed ^ 0x5DEECE66DLL) & ((1LL << 48) - 1);

				int bits = 0, val = 0;
				do {
					seed = (seed * 0x5DEECE66DLL + 0xBLL) & ((1LL << 48) - 1);
					bits = seed >> 17;
					val = bits - bits / 10 * 10;
				} while (bits - val + 9 < 0);

				bool isSlimeChunk = val == 0;
				//-----={End check if slime chunk}=-----//
				if (isSlimeChunk) {
					if (pattern[i * (this->patternWidth - 1) + i] != 1)isNotMatch = true;
					++startSum;
				} else {
					if (pattern[i * (this->patternWidth - 1) + i] != 0)isNotMatch = true;
				}
			}

			if (isNotMatch) {
				Xskip = matchTable[startSum];
				continue;
			}

			bool isDone = false;
			for(int x = this->patternWidth - 2; x != -1; --x) {
				for(int z = 0; z < this->patternHeight; ++z) {
					//-----={Start check if slime chunk}=-----//
					int xr = chunkX + x;
					int zr = chunkZ + z;

					long long seed =
						gameSeed +
						(long long)(xr * xr * 4987142) +
						(long long)(xr * 5947611) +
						((long long)(zr * zr) * (long long)4392871) +
						(long long)(zr * 389711) ^ (long long)987234911;

					seed = (seed ^ 0x5DEECE66DLL) & ((1LL << 48) - 1);

					int bits = 0, val = 0;
					do {
						seed = (seed * 0x5DEECE66DLL + 0xBLL) & ((1LL << 48) - 1);
						bits = seed >> 17;
						val = bits - bits / 10 * 10;
					} while (bits - val + 9 < 0);

					bool isSlimeChunk = val == 0;
					//-----={End check if slime chunk}=-----//
					if(isSlimeChunk != pattern[z * this->patternWidth + x]) {
						Xskip = matchTable[startSum];
						isDone = true;
						break;
					}
				}
				if(isDone)break;
			}
			if(isDone)continue;
			cout << "Found pattern with seed " << gameSeed << " at [" << chunkX * 16 << ", " << chunkZ * 16 << "]." << endl;
			++nmOfPatternsFound;
			if(!searchForAll) {
				auto endTime = chrono::system_clock::now();
				chrono::duration<double> searchDuration = endTime - startTime;
				cout << "Search ended after " << searchDuration.count() << " seconds.\n\n";
				return;
			}
		}
	}

	auto endTime = chrono::system_clock::now();
	chrono::duration<double> searchDuration = endTime - startTime;
	cout << "Search ended after " << searchDuration.count() << " seconds. Found " << nmOfPatternsFound << " matches.\n\n";
}

void PatternSearch::startSearchBedrock(long long gameSeed, vector<pair<int, int>> boundaries, bool searchForAll, bool showDebugInfo) {
	if (showDebugInfo) {
		this->printPattern();
		this->printMatchTable();
		this->printChanceTable();
		this->printAverageSkip();
	}

	int minX = boundaries.at(0).first;
	int maxX = boundaries.at(0).first;
	int minZ = boundaries.at(0).second;
	int maxZ = boundaries.at(0).second;

	int x = boundaries.at(1).first;
	int z = boundaries.at(1).second;
	if (minX > x)minX = x;
	if (maxX < x)maxX = x;
	if (minZ > z)minZ = z;
	if (maxZ < z)maxZ = z;

	int nmOfPatternsFound = 0;
	auto startTime = chrono::system_clock::now();

	int Xskip = 1;
	for (int chunkZ = minZ; chunkZ != maxZ + 1; ++chunkZ) {
		for (int chunkX = minX; chunkX < maxX + 1; chunkX += Xskip) {
			int startSum = 0;
			bool isNotMatch = false;

			for (int i = this->patternHeight - 1; i != -1; --i) {
				//-----={Start check if slime chunk}=-----//
				int x = chunkX + this->patternWidth - 1;
				int z = chunkZ + i;

				int seed = (x * 522133279) ^ z;

				int N = 624;
				int M = 397;
				int* mt = new int[N]; 
				int mti; 
				int mtiFast;
				int UPPER_MASK = 0x80000000;
				int LOWER_MASK = 0x7fffffff;
				int MATRIX_A = 0x9908b0df;
				int MAG_01[] = { 0, MATRIX_A };

				mt[0] = seed;
				for (mtiFast = 1; mtiFast <= M; mtiFast++) {
					mt[mtiFast] = 1812433253
						* ((((unsigned int)mt[mtiFast - 1]) >> 30) ^ mt[mtiFast - 1])
						+ mtiFast;
				}
				mti = N;

				if(mti == N) {
					mti = 0;
				}
				else if(mti > N) {
					mt[0] = 5489;
					for (mtiFast = 1; mtiFast <= M; mtiFast++) {
						mt[mtiFast] = 1812433253
							* (((unsigned int)(mt[mtiFast - 1]) >> 30) ^ mt[mtiFast - 1])
							+ mtiFast;
					}
					mti = 0;
				}

				if (mti >= N - M) {
					if (mti >= N - 1) {
						mt[N - 1] = MAG_01[mt[0] & 1]
							^ (((unsigned int)(mt[0] & LOWER_MASK | mt[N - 1] & UPPER_MASK)) >> 1)
							^ mt[M - 1];
					}
					else {
						mt[mti] = MAG_01[mt[mti + 1] & 1]
							^ (((unsigned int)(mt[mti + 1] & LOWER_MASK | mt[mti] & UPPER_MASK)) >> 1)
							^ mt[mti - (N - M)];
					}
				}
				else {
					mt[mti] = MAG_01[mt[mti + 1] & 1]
						^ (((unsigned int)(mt[mti + 1] & LOWER_MASK | mt[mti] & UPPER_MASK)) >> 1)
						^ mt[mti + M];

					if (mtiFast < N) {
						mt[mtiFast] = 1812433253
							* ((((unsigned int)mt[mtiFast - 1]) >> 30) ^ mt[mtiFast - 1])
							+ mtiFast;
						mtiFast++;
					}
				}

				int ret = mt[mti++];
				ret = ((ret ^ ((unsigned int)ret >> 11)) << 7) & 0x9d2c5680 ^ ret ^ ((unsigned int)ret >> 11);
				ret = (ret << 15) & 0xefc60000 ^ ret ^ (((unsigned int)((ret << 15) & 0xefc60000 ^ ret)) >> 18);

				bool isSlimeChunk = ((unsigned long long)ret) % 10 == 0;
				//-----={End check if slime chunk}=-----//
				if (isSlimeChunk) {
					if (pattern[i * (this->patternWidth - 1) + i] != 1)isNotMatch = true;
					++startSum;
				}
				else {
					if (pattern[i * (this->patternWidth - 1) + i] != 0)isNotMatch = true;
				}
			}

			if (isNotMatch) {
				Xskip = matchTable[startSum];
				continue;
			}

			bool isDone = false;
			for (int x = this->patternWidth + -2; x != -1; --x) {
				for (int z = 0; z < this->patternHeight; ++z) {
					//-----={Start check if slime chunk}=-----//
					int xr = chunkX + x;
					int zr = chunkZ + z;

					int seed = (xr * 522133279) ^ zr;

					int N = 624;
					int M = 397;
					int* mt = new int[N];
					int mti;
					int mtiFast;
					int UPPER_MASK = 0x80000000;
					int LOWER_MASK = 0x7fffffff;
					int MATRIX_A = 0x9908b0df;
					int MAG_01[] = { 0, MATRIX_A };

					mt[0] = seed;
					for (mtiFast = 1; mtiFast <= M; mtiFast++) {
						mt[mtiFast] = 1812433253
							* ((((unsigned int)mt[mtiFast - 1]) >> 30) ^ mt[mtiFast - 1])
							+ mtiFast;
					}
					mti = N;

					if (mti == N) {
						mti = 0;
					}
					else if (mti > N) {
						mt[0] = 5489;
						for (mtiFast = 1; mtiFast <= M; mtiFast++) {
							mt[mtiFast] = 1812433253
								* (((unsigned int)(mt[mtiFast - 1]) >> 30) ^ mt[mtiFast - 1])
								+ mtiFast;
						}
						mti = 0;
					}

					if (mti >= N - M) {
						if (mti >= N - 1) {
							mt[N - 1] = MAG_01[mt[0] & 1]
								^ (((unsigned int)(mt[0] & LOWER_MASK | mt[N - 1] & UPPER_MASK)) >> 1)
								^ mt[M - 1];
						}
						else {
							mt[mti] = MAG_01[mt[mti + 1] & 1]
								^ (((unsigned int)(mt[mti + 1] & LOWER_MASK | mt[mti] & UPPER_MASK)) >> 1)
								^ mt[mti - (N - M)];
						}
					}
					else {
						mt[mti] = MAG_01[mt[mti + 1] & 1]
							^ (((unsigned int)(mt[mti + 1] & LOWER_MASK | mt[mti] & UPPER_MASK)) >> 1)
							^ mt[mti + M];

						if (mtiFast < N) {
							mt[mtiFast] = 1812433253
								* ((((unsigned int)mt[mtiFast - 1]) >> 30) ^ mt[mtiFast - 1])
								+ mtiFast;
							mtiFast++;
						}
					}

					int ret = mt[mti++];
					ret = ((ret ^ ((unsigned int)ret >> 11)) << 7) & 0x9d2c5680 ^ ret ^ ((unsigned int)ret >> 11);
					ret = (ret << 15) & 0xefc60000 ^ ret ^ (((unsigned int)((ret << 15) & 0xefc60000 ^ ret)) >> 18);

					bool isSlimeChunk = ((unsigned long long)ret) % 10 == 0;
					//-----={End check if slime chunk}=-----//
					if (isSlimeChunk != pattern[z * this->patternWidth + x]) {
						Xskip = matchTable[startSum];
						isDone = true;
						break;
					}
				}
				if (isDone)break;
			}
			if (isDone)continue;
			cout << "Found pattern with seed " << gameSeed << " at [" << chunkX * 16 << ", " << chunkZ * 16 << "]." << endl;
			++nmOfPatternsFound;
			if (!searchForAll) {
				auto endTime = chrono::system_clock::now();
				chrono::duration<double> searchDuration = endTime - startTime;
				cout << "Search ended after " << searchDuration.count() << " seconds.\n\n";
				return;
			}
		}
	}

	auto endTime = chrono::system_clock::now();
	chrono::duration<double> searchDuration = endTime - startTime;
	cout << "Search ended after " << searchDuration.count() << " seconds. Found " << nmOfPatternsFound << " matches.\n\n";
}

//--------------={PRECOMPUTATION STARTS HERE}=--------------//

void PatternSearch::preComputeMatchTable() {
	matchTable = new int[this->patternHeight + 1];

	for(int i = 0; i != this->patternHeight + 1; ++i) {
		matchTable[i] = this->patternWidth;
	}

	for(int w = 0; w != this->patternWidth; ++w) {
		int sum = 0;
		for(int h = 0; h != this->patternHeight; ++h) {
			sum += pattern[h * this->patternWidth + w] == 1;
		}
		matchTable[sum] = this->math->getMax(1, this->patternWidth - w - 1);
	}
}

void PatternSearch::preComputeChanceTable() {
	chanceTable = new double[this->patternHeight + 1];

	for(int t = 0; t < this->patternHeight + 1; ++t) {
		double sum = 1;
		for(int i = 0; i < t; ++i)sum *= 0.1;
		for(int j = 0; j < this->patternHeight - t; ++j)sum *= 0.9;

		double multiplier = this->math->getFactorial(this->patternHeight) / 
							(this->math->getFactorial(t) * this->math->getFactorial(this->patternHeight - t));
		sum *= multiplier;
		chanceTable[t] = sum;
	}
}

void PatternSearch::preComputeAverageSkip() {
	float average = 0;
	for(int i = 0; i != this->patternHeight + 1; ++i) {
		average += this->matchTable[i] * this->chanceTable[i];
	}
	this->averageSkip = average;
}

//--------------={DEBUG PRINTING STUFF STARTS HERE}=--------------//

void PatternSearch::printPattern() {
	cout << "The pattern for this search is as follows.\n\n";

	for(int h = 0; h != this->patternHeight; ++h) {
		for(int w = 0; w != this->patternWidth; ++w) {
			cout << (this->pattern[h * this->patternWidth + w] == 1 ? "1" : "0");
		}
		cout << "\n";
	}

	cout << "\n";
}

void PatternSearch::printMatchTable() {
	for(int i = 0; i != this->patternHeight + 1; ++i) {
		cout << "Mistmatch at a column with " << i << (i == 1 ? " chunk " : " chunks ") << "calls for ";
		cout << matchTable[i] << (matchTable[i] == 1 ? " skip." : " skips.") << "\n";
	}

	cout << "\n";
}

void PatternSearch::printChanceTable() {
	for(int i = 0; i < this->patternHeight + 1; ++i) {
		cout << "The chance of finding " << i << (i == 1 ? " slime chunk" : " slime chunks") << " in a given " << this->patternHeight << " chunks is ";
		cout << this->chanceTable[i] * 100 << "%";
		long long possibilities = this->math->getFactorial(this->patternHeight) / 
									(long long)(this->math->getFactorial(i) * this->math->getFactorial(this->patternHeight - i));
		cout << " accross " << possibilities << (possibilities == 1 ? " possiblity." : " possibilities.");
		cout << "\n";
	}

	cout << "\n";
}

void PatternSearch::printAverageSkip() {
	cout << "\n";
	cout << "In average, there is " << 0.1 * this->patternHeight << " slime chunks in a given " << this->patternHeight << " chunks.\n";
	cout << "Therefore, the average skip per scan is " << this->averageSkip << ", with a maximum of ";

	int max = 0;
	int min = this->patternHeight + 1;

	for(int i = 0; i < this->patternHeight + 1; ++i) {
		if(this->matchTable[i] > max)max = this->matchTable[i];
		if(this->matchTable[i] < min)min = this->matchTable[i];
	}

	cout << max << " and a minimum of " << min << ".\n\n";
}
