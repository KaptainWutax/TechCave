#include "Math.h"
#include <iostream>

Math::Math(int factorialPrecision) {
	preComputeFactorialTable(factorialPrecision);
}

Math::~Math() {
}

void Math::preComputeFactorialTable(int factorialPrecision) {
	factorialTable = new double[factorialPrecision];

	for(int f = 0; f != factorialPrecision + 1; ++f) {
		factorialTable[f] = factorial(f);
	}
}

double Math::factorial(int n) {
	return (n > 1 ? n * factorial(n - 1) : 1);
}
