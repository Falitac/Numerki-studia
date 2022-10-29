#include <iostream>
#include <vector>
#include <fstream>

using Matrix = std::vector<std::vector<float>>;

void gauss(Matrix& a, std::vector<float>& b, std::vector<float>&x) {
	int n = a.size();
	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			float m = a[j][i] / a[i][i];
			b[j] -= b[i] * m;
			for (int k = j; k < n; k++) {
				a[j][k] -= a[i][k] * m;
			}
		}
	}
	// obliczanie zmiennych:
	for (int i = n - 1; i >= 0; i--) {
		x[i] = b[i];
		for (int j = n - 1; j > i; j--) {
			x[i] -= a[i][j] * x[j];
		}
		x[i] /= a[i][i];
	}
}

void print(std::vector<float>& a) {
	for (int i = 0; i < a.size(); i++) {
		std::printf("[%d]: %5.2f\n", i, a[i]);
	}
}

void print(Matrix& a) {
	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < a[i].size(); j++) {
			std::printf("%5.2f ", a[i][j]);
		}
		std::printf("\n");
	}
}

void jacob(Matrix& a, std::vector<float>& b, std::vector<float>& x, int iters) {
	int n = a.size();
	for (int it = 0; it < iters; it++) {

	}
}

bool isWeakDominant(Matrix& a) {
	int n = a.size();
	bool good = false;
	for (int i = 0; i < n; i++) {
		float sum = 0.f;
		for (int j = 0; j < n; j++) {
			if (i == 0) continue;
			sum += std::fabs(a[i][j]);
		}
		if (sum > std::fabs(a[i][i])) {
			return false;
		}
		if (sum > std::fabs(a[i][i])) {
			good = true;
		}
	}
	return good;
}

int main() {
	std::fstream input("RURL_dane2.txt");
	if (!input.good()) {
		std::printf("nope\n");
		return 1;
	}
	int n;
	input >> n;
	Matrix m = std::vector<std::vector<float>>(n, std::vector<float>(n, 0.0f));

	std::vector<float> b(n);
	std::vector<float> x(n, 0.f);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			input >> m[i][j];
			if (i == j && std::fabs(m[i][j]) < FLT_EPSILON) {
				std::printf("Error matrix 0. at [%d][%d]\n", i, i);
				return 1;
			}
		}
		input >> b[i];
	}

	std::printf("Macierz [%i][%i]:\n", n, n);
	print(m);
	std::printf("Wektor wyr. wolnych:\n");
	print(b);
	std::printf("Czy s³abo dom.?: %s", isWeakDominant(m) ? "Tak" : "Nie");

	std::puts("Po obl.:\n");
	print(m);
	std::printf("Wektor wyr. wolnych:\n");
	print(b);
	std::printf("Wynik:\n");
	print(x);

}