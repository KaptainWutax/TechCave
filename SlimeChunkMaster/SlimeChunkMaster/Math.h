#pragma once
class Math {
private:
	double* factorialTable;

	void preComputeFactorialTable(int factorialPrecision);
	double factorial(int n);
public:
	Math(int factorialPrecision);
	~Math();

	template<class T, class S>
	T getMax(T a, S b);

	template<class T>
	double getFactorial(T f);
};

//A template method definition has to be in the header...what? LOL xD
template<class T, class S>
T Math::getMax(T a, S b) {
	return (a > b ? a : b);
}

template<class T>
double Math::getFactorial(T f) {
	return factorialTable[(int)f];
}
