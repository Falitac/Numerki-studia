#include <iostream>
#include <vector>
#include <fstream>

using Matrix = std::vector<std::vector<float>>;

// Klasa pomocnicza do reprezentacji macierzy. 
namespace fil {
	class Matrix {
	public:
		Matrix(size_t row, size_t col, float val = 0.f)
		: data(row * col, val)
		, row(row)
		, col(col) {
		}
		float* operator[](size_t n) {
			return &data[n * col];
		}
	private:
		size_t row;
		size_t col;
		std::vector<float> data;
	};
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
void gaussEliminateToUpper(Matrix& a, std::vector<float>& b) {
	int n = a.size();
	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			float m = a[j][i] / a[i][i];
			b[j] -= b[i] * m;
			for (int k = i; k < n; k++) {
				a[j][k] -= a[i][k] * m;
			}
		}
	}
}

// obliczanie zmiennych:
std::vector<float> solveUpper(Matrix& a, std::vector<float>& b) {
	int n = a.size();
	std::vector<float> x(n);

	for (int i = n - 1; i >= 0; i--) {
		x[i] = b[i];
		for (int j = n - 1; j > i; j--) {
			x[i] -= a[i][j] * x[j];
		}
		x[i] /= a[i][i];
	}
	return x;
}

std::vector<float> jacob(Matrix& a, std::vector<float>& b, int iters) {
	int n = a.size();
	Matrix LU = a;
	std::vector<float> x(n, 0.f);
	std::vector<float> d(n);
	for (int i = 0; i < n; i++) {
		d[i] = 1 / a[i][i];
		LU[i][i] = 0.f;
	}
	std::printf("D:\n");
	print(d);
	std::printf("LU:\n");
	print(LU);

	Matrix dlu = std::vector<std::vector<float>>(n, std::vector<float>(n, 0.0f));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			dlu[i][j] = -d[i] * LU[i][j];
		}
	}
	std::printf("DLU:\n");
	print(dlu);
	std::vector<float> db = b;
	for (int i = 0; i < n; i++) {
		db[i] *= d[i];
	}
	std::printf("D^-1*B:\n");
	print(db);
	for (int it = 0; it < iters; it++) {
		//x_i+1 = -(d * LU) * x_i + d*b
		std::vector<float> tempx = x;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				x[i] += dlu[i][j] * tempx[i];
			}
			x[i] += db[i];
		}
	}
	return x;
}

int main() {
	std::fstream input("test2.txt");
	if (!input.good()) {
		std::printf("nope\n");
		return 1;
	}
	int n;
	input >> n;
	Matrix m = std::vector<std::vector<float>>(n, std::vector<float>(n, 0.0f));

	std::vector<float> b(n);
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
	std::printf("Czy slabo dom.?: %s\n", isWeakDominant(m) ? "Tak" : "Nie");

	//gaussEliminateToUpper(m, b);
	//auto x = solveUpper(m, b);
	auto x = jacob(m, b, 5);

	std::puts("Macierz U:");
	print(m);
	std::printf("Wektor wyr. wolnych:\n");
	print(b);
	std::printf("Wynik:\n");
	print(x);
}