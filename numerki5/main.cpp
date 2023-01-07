#include <iostream>
#include <vector>
#include <format>

struct Node {
	double x;
	double weight = 1.;
};

double integrateRectMethod(auto foo, double a, double b, int n) {
	double result = 0.0;
	double width = (b - a) / static_cast<double>(n);
	for (int i = 0; i < n; i++) {
		double x = width * (i + 0.5) + a;
		result += foo(x);
	}
	result *= width;

	return result;
}

double integrateTrapMethod(auto foo, double a, double b, int n) {
	double result = 0.0;
	double width = (b - a) / static_cast<double>(n);
	for (int i = 0; i < n; i++) {
		double x = width * i + a;
		result += (foo(x) + foo(x + width));
	}
	result *= width / 2.;

	return result;
}

double integrateSimpMethod(auto foo, double a, double b, int n) {
	double result = 0.0;
	double width = (b - a) / static_cast<double>(n);
	for (int i = 0; i < n; i++) {
		double x = width * i + a;
		result += (foo(x) + 4. * foo(x + width * .5) + foo(x + width));
	}
	result *= width / 6.;

	return result;
}

double integrateGaLageMethod(auto foo, double a, double b, std::vector<Node>& nodes) {
	int n = nodes.size();
	double result = 0.0;
	double width = (b - a) / static_cast<double>(n);

	for (int i = 0; i < n; i++) {
		double x = (a + b) * .5 + (b - a) * .5 * nodes[i].x;
		result += nodes[i].weight * foo(x);
	}
	result *= (b - a) * .5;

	return result;
}

void proceedFunction(auto foo, double a, double b, int iters = 5) {
	auto res1 = integrateRectMethod(foo, a, b, iters);
	auto res2 = integrateTrapMethod(foo, a, b, iters);
	auto res3 = integrateSimpMethod(foo, a, b, iters);

	std::printf("Metoda prostakatow: %.8f\n", res1);
	std::printf("Metoda trapezow:    %.8f\n", res2);
	std::printf("Metoda Simpsona:    %.8f\n", res3);
}

int main(int argc, char** argv) {
	auto foo1 = [](double x) {
		return std::sin(x);
	};
	auto foo2 = [](double x) {
		return x*x + 2*x + 5;
	};
	auto foo3 = [](double x) {
		return std::exp(x);
	};

	std::vector<Node> nodes0 = {
		{-0.57735, 1.},
		{ 0.57735, 1.},
	};
	std::vector<Node> nodes1 = {
		{-0.861136, 0.347855},
		{-0.339981, 0.652145},
		{ 0.339981, 0.652145},
		{ 0.861136, 0.347855},
	};

	std::cout << "f(x) = sin(x), 0.5 <= x <= 2.5\n";
	proceedFunction(foo1, .5, 2.5);

	auto resGaLage1 = integrateGaLageMethod(foo1, .5, 2.5, nodes0);
	auto resGaLage2 = integrateGaLageMethod(foo1, .5, 2.5, nodes1);

	printf("Metoda Gaussa-Legendre'a (2 wezly): %.8f\n", resGaLage1);
	printf("Metoda Gaussa-Legendre'a (4 wezly): %.8f\n", resGaLage2);

	std::cout << "f(x) = x^2 + 2x + 5, 0.5 <= x <= 5\n";
	proceedFunction(foo2, .5, 5);

	std::cout << "f(x) = exp(x), 0.5 <= x <= 5\n";
	proceedFunction(foo3, .5, 5);


	return 0;
}