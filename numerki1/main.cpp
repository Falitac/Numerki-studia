#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cfloat>


// Klasa pomocnicza do reprezentacji macierzy. 
class Matrix {
public:
	Matrix(size_t row, size_t col, float val = 0.f)
	: row(row)
	, col(col)
	, data(row * col, val) {
	}

	float* operator[](size_t n) {
		return &data[n * col];
	}
	size_t rowSize() const { return row; };
	size_t colSize() const { return col; };
private:
	size_t row;
	size_t col;
	std::vector<float> data;
};

void print(std::vector<float>& a) {
	for (int i = 0; i < a.size(); i++) {
		std::printf("[%d]: %6.5g\n", i, a[i]);
	}
}

void print(Matrix& a) {
	for (int i = 0; i < a.rowSize(); i++) {
		for (int j = 0; j < a.colSize(); j++) {
			std::printf("%6.5g ", a[i][j]);
		}
		std::printf("\n");
	}
}

bool isWeakDominant(Matrix& a) {
	int n = a.rowSize();
	bool good = false;
	for (int i = 0; i < n; i++) {
		float sum = 0.f;
		for (int j = 0; j < n; j++) {
			if (i == j) continue;
			sum += std::fabs(a[i][j]);
		}

		if (std::fabs(a[i][i]) < sum) {
			return false;
		} else {
			good = true;
		}
	}
	return good;
}

int pivotOnDiagonal(Matrix& a, std::vector<float>& b, int diagonal) {
  int highest = diagonal;
  float maxVal = std::fabs(a[diagonal][diagonal]);
  for(int i = diagonal + 1; i < a.rowSize(); i++) {
    auto absVal = std::fabs(a[i][diagonal]);
    if(absVal > maxVal) {
      maxVal = absVal;
      highest = i;
    }
  }
  return highest;
}

void swapRows(Matrix& a, std::vector<float>& b, int rowIndex1, int rowIndex2) {
  if(rowIndex1 == rowIndex2) return;
  std::swap(b[rowIndex1], b[rowIndex2]);
  std::swap_ranges(a[rowIndex1], a[rowIndex1] + a.colSize(), a[rowIndex2]);
}

void gaussEliminateToUpper(Matrix& a, std::vector<float>& b) {
	int n = a.rowSize();
	for (int i = 0; i < n - 1; i++) {
    int highest = pivotOnDiagonal(a, b, i);
    swapRows(a, b, i, highest);
    if(std::fabs(a[i][i]) <= 1e-5) {
      std::printf("well, 0");
      throw std::exception();
    }
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
	int n = a.rowSize();
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

std::vector<float> solveByJacobi(Matrix a, std::vector<float> b, int iters) {
	int n = a.rowSize();
	std::vector<float> x(n, 0.f);

	std::vector<float> d(n);
	for (int i = 0; i < n; i++) {
		d[i] = 1 / a[i][i];
		a[i][i] = 0.f;
	}
  
  std::puts("Macierz L+U:");
  print(a);
  std::puts("Macierz odwrotna diagonalna D^(-1):\n(reprezentacja w wektorze)");
  print(d);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			a[i][j] = -d[i] * a[i][j];
		}
		b[i] *= d[i];
	}
	for (int it = 0; it < iters; it++) {
		std::vector<float> tempx(n, 0.f);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				tempx[i] += a[i][j] * x[j];
			}
			tempx[i] += b[i];
		}
    x = tempx;
	}
	return x;
}

std::pair<Matrix, std::vector<float>> loadMatrixWithVector(const std::string& fileName) {
	std::fstream input(fileName);
	if (!input.good()) {
		std::printf("Nope, problem z plikiem\n");
    std::exit(1);
	}

	int n;
	input >> n;
	Matrix m(n, n);
	std::vector<float> b(n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			input >> m[i][j];
			if (i == j && std::fabs(m[i][j]) < FLT_EPSILON) {
				std::printf("Error matrix 0. at [%d][%d]\n", i, i);
        std::exit(1);
			}
		}
		input >> b[i];
	}
  return {m, b};
}

void gauss() {
  auto [m, b] = loadMatrixWithVector("test4.txt");
	std::printf("Macierz %ix%i:\n", m.rowSize(), m.rowSize());
	print(m);
	std::printf("Wektor wyr. wolnych:\n");
	print(b);

  try {
    gaussEliminateToUpper(m, b);
    print(m);
  } catch(std::exception& e) {
    std::printf("%s\n", e.what());
  }
  auto x = solveUpper(m, b);
	std::printf("Wynik:\n");
	print(x);
}

void jacobi() {
  auto [m, b] = loadMatrixWithVector("test4.txt");

	std::printf("Macierz %ix%i:\n", m.rowSize(), m.rowSize());
	print(m);
	std::printf("Wektor wyr. wolnych:\n");
	print(b);
	std::printf("Czy slabo dom.?: %s\n", isWeakDominant(m) ? "Tak" : "Nie");

	//gaussEliminateToUpper(m, b);
	//auto x = solveUpper(m, b);
	int iters = 1;
	puts("Podaj ilość iteracji:\n");
	std::cin >> iters;
	auto x = solveByJacobi(m, b, iters);
	std::printf("Wynik:\n");
	print(x);
}

int main() {
  gauss();
  jacobi();
}

