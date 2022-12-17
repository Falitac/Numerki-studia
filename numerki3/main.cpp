#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

struct Point {
	double x;
	double y;
};

std::vector<Point> loadPoints(const std::string& path) {
	std::fstream input(path);
	int n;
	input >> n;

	std::vector<Point> result(n);
	for (int i = 0; i < n; i++) {
		input >> result[i].x >> result[i].y;
	}

	return result;
}

double langrange(std::vector<Point>& points, double x) {
	double result = 0.0;

	for (int i = 0; i < points.size(); i++) {
		double res = 1.0;
		for (int j = 0; j < points.size(); j++) {
			if (i == j) {
				continue;
			}
			res *= (x - points[j].x) / (points[i].x - points[j].x);
		}
		res *= points[i].y;
		result += res;
	}
	return result;
}

void printVec(std::vector<double>& vec, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%f ", vec[n * i + j]);
		}
		puts("");
	}
}

double newton(std::vector<Point>& points, double x) {
	double result = points[0].y;
	double p = 1.0;
	std::printf("b0: %f\n", result);
	for (int i = 1; i < points.size(); i++) {
		double b = 0.;
		for (int j = 0; j <= i; j++) {
			double den = points[j].y;
			for (int k = 0; k <= i; k++) {
				if (k == j) {
					continue;
				}
				den /= (points[j].x - points[k].x);
			}
			b += den;
		}
		p *= x - points[i - 1].x;
		std::printf("b%i: %f\n", i, b);
		result += b * p;
	}
	return result;
}

double newton2(std::vector<Point>& points, double x) {
	double result = 0.0;
	int n = points.size();
	std::vector<double> diffs(n * n, 0.);

	for (int i = 0; i < n; i++) {
		diffs[i * n] = points[i].y;
	}

	for (int i = 1; i < n; i++) {
		for (int j = 0; j < n - i; j++) {
			diffs[j * n + i] = diffs[(j + 1) * n + i - 1] - diffs[(j)*n + i - 1];
			diffs[j * n + i] /= points[j + i].x - points[j].x;
		}
	}

	for (int i = 0; i < n; i++) {
		std::printf("b%d: %7.3f\n", i, diffs[i]);
	}

	double p = 1.0;
	for (int i = 0; i < n; i++) {
		result += p * diffs[i];
		p *= x - points[i].x;
	}
	return result;
}


void printPointsInfo(const std::vector<Point>& points) {
	std::cout << "Liczba wezlow:" << points.size() << '\n';
	for (int i = 0; i < points.size(); i++) {
		std::printf("%i: %.4f %.4f\n", i, points[i].x, points[i].y);
	}
}

int main(int argc, char** argv) {
	auto points = loadPoints("data1.txt");
	printPointsInfo(points);

	std::cout << "Podaj punkt (argument): ";
	double x;
	std::cin >> x;

	std::cout << "Metoda? (L, N, domyślnie langrange'a)\n";
	char choose = 'L';
	std::cin >> choose;

	if (choose == 'N') {
		double y = newton(points, x);
		std::printf("W(%lf) = %lf", x, y);
		return 0;
	}
	double y = langrange(points, x);
	std::printf("L(%lf) = %lf", x, y);

	return 0;
}

