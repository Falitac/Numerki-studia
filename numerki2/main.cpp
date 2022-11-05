#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <cmath>
#include <cfloat>

constexpr double EPS = .01;


template<typename F>
double bisection(F foo, double a, double b) {
	std::puts("i, a, b, foo(a), foo(b), x0, foo(x0), foo(a) * foo(b)");
	for (int i = 0;; i++) {
		double x0 = (b - a) / 2. + a;
		std::printf("%i %7.7lf %7.7lf %7.7lf %7.7lf %7.7lf %7.7lf %7.7lf\n",
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
	for (int i = 0;; i++) {
		double x0 = a - foo(a) * (b - a) / (foo(b) - foo(a));

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
double newton(F foo, Fd der, double x0) {
	std::puts("i x0,  f(x0), f(x0) / f\'(x0)");
	for (int i = 0;; i++) {
		x0 = x0 - foo(x0) / der(x0);
		std::printf("%i %7.7lf %7.7lf %7.7lf\n",
			i, x0, foo(x0), foo(x0) / der(x0));
		if (std::fabs(foo(x0)) < EPS) {
			return x0;
		}
	}
}

template<typename F>
double secant(F foo, double x1, int iters) {
	std::printf("%7s %7s\n", "x0",  "f(x0)");
	double x0 = x1 - .1;
	for (int i = 0; i < iters; i++) {
		std::printf("%i ", i);
		std::printf("%7.7lf %7.7lf %7.7lf\n", x0, x1, foo(x0));
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
		return 2.4 - 1.72 * x + 0.018 * x * x + 0.096 * x * x * x - 0.01 * std::pow(x, 4);
	};
	std::puts("Wlasna funkcja:");
	std::puts("f(x) = -x*(x-8)*(x-5)*(x+6)*(x-5) * .002");
	double x0;
	std::puts("Podaj x0:");
	std::cin >> x0;

	int iters;
	std::puts("Podaj ilosc iteracji:");
	std::cin >> iters;


	std::printf("Bisection: %.3f\n", bisection(ownFun, 6.5, 10.));
	std::printf("Falsi: %.3f\n", falsi(ownFun, 6.5, 10.));
	double newtonResult = newton(ownFun, ownFunD, x0);
	double secantResult = secant(ownFun, x0, iters);
	std::printf("Newton: %.3f\n", newtonResult);
	std::printf("Secant: %.3f\n", secantResult);
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
	std::printf("Bisection: %.3f\n", bisection(foo1, 0., 2.));
	std::printf("Falsi: %.3f\n", falsi(foo1, 0., 2.));


	std::puts("f(x) = -x^3 + 10x + 5.");

	double x0;
	std::puts("Podaj x0:");
	std::cin >> x0;

	double newtonResult = newton(foo2, foo2D, x0);
	std::printf("Newton: %.3f\n", newtonResult);

	int iters;
	std::puts("Podaj ilosc iteracji:");
	std::cin >> iters;

	double secantResult = secant(foo2, x0, iters);
	std::printf("Secant: %.3f\n", secantResult);
	std::printf("Roznica pomiedzy metodami: %.20lf\n", std::fabs(secantResult - newtonResult));

	ownTest();

}