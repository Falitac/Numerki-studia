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
	std::vector<double> quotients(points.size() * points.size(), 0.0);
	double result = 0.0;
	for (int i = 0; i < quotients.size(); i++) {
		quotients[i] = points[i].y;
	}
	for (int i = 1; i < points.size(); i++) {
		for (int j = 0; j < points.size() - i; j++) {
			double up = (quotients[points.size() * (i - 1) + j + 1] - quotients[points.size() * (i - 1) + j]);
			double down = points[j + 1].x - points[j].x;
			quotients[points.size() * i + j] = up / down;
		}
	}
	printVec(quotients, points.size());

	return result;
}

void printPointsInfo(const std::vector<Point>& points) {
	std::cout << "Liczba wezlow:" << points.size() << '\n';
	for (int i = 0; i < points.size(); i++) {
		std::printf("%i: %.4f %.4f\n", i, points[i].x, points[i].y);
	}
}

int main1(int argc, char** argv) {
	auto points = loadPoints("data0.txt");
	printPointsInfo(points);

	std::cout << "Podaj punkt (argument): ";
	double x;
	std::cin >> x;
	double y = langrange(points, x);

	std::printf("L(%lf) = %lf", x, y);

	return 0;
}

int main(int argc, char** argv) {
	auto points = loadPoints("data1.txt");
	printPointsInfo(points);

	std::cout << "Podaj punkt (argument): ";
	double x;
	std::cin >> x;
	double y = newton(points, x);

	std::printf("W(%lf) = %lf", x, y);

	return 0;
}