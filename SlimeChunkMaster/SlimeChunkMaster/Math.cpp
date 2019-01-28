#include "Math.h"

Math::Math(int factorialPrecision) {
	this->preComputeFactorialTable(factorialPrecision);
}

Math::~Math() {
}

void Math::preComputeFactorialTable(int factorialPrecision) {
	this->factorialTable = new double[factorialPrecision];

	for(int f = 0; f != factorialPrecision + 1; ++f) {
		this->factorialTable[f] = this->factorial(f);
	}
}

double Math::factorial(int n) {
	return (n > 1 ? n * this->factorial(n - 1) : 1);
}
