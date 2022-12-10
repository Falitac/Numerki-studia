#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <cmath>
#include <cfloat>

constexpr double EPS = .01;


template<typename F>
double bisection(F foo, double a, double b) {
	std::puts("i\ta\t\tb\t\tf(a)\t\tf(b)\t\tx0\t\tf(x0)\t\tf(a) * f(b)");
	for (int i = 0;; i++) {
		double x0 = (b - a) / 2. + a;
		std::printf("%i\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
			          i, a, b, foo(a), foo(b), x0, foo(x0), foo(a) * foo(b));
		if (std::fabs(foo(x0)) < EPS) {
			return x0;
		}
		if (foo(a) * foo(x0) < 0) {
			b = x0;
		} else {
			a = x0;
		}
	}
	return 0.;
}

template<typename F>
double falsi(F foo, double a, double b) {
	std::puts("i\ta\t\tb\t\tx0\t\tf(x)");
	for (int i = 0;; i++) {
		double x0 = a - foo(a) * (b - a) / (foo(b) - foo(a));
		std::printf("%i\t%7.8lf\t%7.8lf\t%7.8lf\t%7.8f\n",
			i, a, b, x0, foo(x0));

		if (std::fabs(foo(x0)) < EPS) {
			return x0;
		}
		if (foo(a) * foo(x0) < 0) {
			b = x0;
		}
		else {
			a = x0;
		}
	}
}

template<typename F, typename Fd>
double newton(F foo, Fd der, double x0, int iters) {
	std::puts("i\tx0,\t\tf(x0),\t\tf(x0)/f'(x0)");
	for (int i = 0; i < iters; i++) {
		x0 = x0 - foo(x0) / der(x0);
		std::printf("%i\t%7.8lf\t%7.8lf\t%7.8lf\n",
			i, x0, foo(x0), foo(x0) / der(x0));
	}
	return x0;
}

template<typename F>
double secant(F foo, double x1, int iters) {
	std::puts("i\tx0\t\tx1\t\tf(x0)");
	double x0 = x1 - .1;
	for (int i = 0; i < iters; i++) {
		std::printf("%i\t%7.8lf\t%7.8lf\t%7.8lf\n",
			i, x0, x1, foo(x1));
		double tmp = x1;
		x1 = x1 - foo(x1) * (x1 - x0) / (foo(x1) - foo(x0));
		x0 = tmp;
	}
	return x0;
}

void ownTest() {
	auto ownFun = [](double x) {
		return -x * (x - 8) * (x - 5) * (x + 6) * (x - 5.5) * .002;
	};
	auto ownFunD = [](double x) {
		auto x2 = x * x;
		return -0.01 * x2 * x2 + 0.1 * x2 * x - 0.003 * x2 - 1.796 * x + 2.64;
	};
	std::puts("Wlasna funkcja:");
	std::puts("f(x) = -x*(x-8)*(x-5)*(x+6)*(x-5.5) * .002");
	double x0;
	std::puts("Podaj x0:");
	std::cin >> x0;

	int iters;
	std::puts("Podaj ilosc iteracji:");
	std::cin >> iters;


	std::printf("Bisection: %.9f\n", bisection(ownFun, 6.5, 10.));
	std::printf("Falsi: %.9f\n", falsi(ownFun, 6.5, 10.));

	double newtonResult = newton(ownFun, ownFunD, x0, iters);
	std::printf("Newton: %.9f\n", newtonResult);
	double secantResult = secant(ownFun, x0, iters);
	std::printf("Secant: %.9f\n", secantResult);
	std::printf("Roznica pomiedzy metodami: %.20lf\n", std::fabs(secantResult - newtonResult));
}

int main() {
	auto foo1 = [](double x) {
		return cos(x*x*x - 2 * x);
	};
	auto foo2 = [](double x) {
		return -x * x * x + 10. * x + 5.;
	};
	auto foo2D = [](double x) {
		return -3 * x * x + 10.;
	};

	std::puts("f(x) = cos(x^3-2x)");

	auto bisectionResult = bisection(foo1, 0., 2.);
	std::printf("Bisection: %.9f\n", bisectionResult);

	auto falsiResult = falsi(foo1, 0., 2.);
	std::printf("Falsi: %.9f\n", falsiResult);
	std::printf("|Falsi - Bisection| = %.9f\n", std::fabs(bisectionResult - falsiResult));


	std::puts("f(x) = -x^3 + 10x + 5");
	std::puts("f'(x) = -3x^2 + 10");

	double x0;
	std::puts("Podaj x0:");
	std::cin >> x0;

	int iters;
	std::puts("Podaj ilosc iteracji:");
	std::cin >> iters;

	double newtonResult = newton(foo2, foo2D, x0, iters);
	std::printf("Newton: %.9f\n", newtonResult);

	double secantResult = secant(foo2, x0, iters);
	std::printf("Secant: %.9f\n", secantResult);
	std::printf("Roznica pomiedzy metodami: %.20lf\n", std::fabs(secantResult - newtonResult));

	ownTest();

}