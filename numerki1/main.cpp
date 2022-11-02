#include <iostream>
#include <algorithm>
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

    void print() {
        for (int i = 0; i < rowSize(); i++) {
            for (int j = 0; j < colSize(); j++) {
                std::printf("%9.8g ", (*this)[i][j]);
            }
            std::printf("\n");
        }
    }
private:
    size_t row;
    size_t col;
    std::vector<float> data;
};

void print(std::vector<float>& a) {
    for (int i = 0; i < a.size(); i++) {
        std::printf("[%d]: %9.8g\n", i, a[i]);
    }
}

bool isWeakDominant(Matrix& a) {
    int n = a.rowSize();
    for (int i = 0; i < n; i++) {
        float sum = 0.f;
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            sum += std::fabs(a[i][j]);
        }

        if (std::fabs(a[i][i]) < sum) {
            return false;
        }
    }
    return true;
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
        for (int j = i + 1; j < n; j++) {
            float m = a[j][i] / a[i][i];
            b[j] -= b[i] * m;
            for (int k = i; k < n; k++) {
                a[j][k] -= a[i][k] * m;
            }
        }
    }
}

void gaussEliminateToUpperPivoting(Matrix& a, std::vector<float>& b) {
    int n = a.rowSize();
    for (int i = 0; i < n - 1; i++) {
        int highest = pivotOnDiagonal(a, b, i);
        swapRows(a, b, i, highest);
        for (int j = i + 1; j < n; j++) {
            float m = a[j][i] / a[i][i];
            b[j] -= b[i] * m;
            for (int k = i; k < n; k++) {
                a[j][k] -= a[i][k] * m;
            }
        }
    }
}

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
    a.print();
    std::puts("Macierz odwrotna diagonalna D^(-1):\n(reprezentacja w wektorze)");
    print(d);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = -d[i] * a[i][j];
        }
        b[i] *= d[i];
    }
    std::vector<float> tempx(n);
    for (int it = 0; it < iters; it++) {
        std::fill(tempx.begin(), tempx.end(), 0.f);
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
            if (i == j && std::fabs(m[i][j]) == 0.f) {
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
    std::printf("Macierz %dx%d:\n", m.rowSize(), m.rowSize());
    m.print();
    std::printf("Wektor wyr. wolnych:\n");
    print(b);

    int choose;
    std::puts("Użyć pivoting'u? (0 - nie, inna wartość - tak)");
    std::cin >> choose;
    if (choose) {
        gaussEliminateToUpperPivoting(m, b);
    } else {
        gaussEliminateToUpper(m, b);
    }
    auto x = solveUpper(m, b);

    std::printf("Po przekształceniu:\n");
    m.print();
    std::printf("Wektor b:\n");
    print(b);
    std::printf("Wynik:\n");
    print(x);
}

void jacobi() {
    auto [m, b] = loadMatrixWithVector("test2.txt");

    std::printf("Macierz %dx%d:\n", m.rowSize(), m.rowSize());
    m.print();
    std::printf("Wektor wyr. wolnych:\n");
    print(b);
    std::printf("Czy slabo dominujaca na diagonali?: %s\n", isWeakDominant(m) ? "Tak" : "Nie");

    int iters = 1;
    puts("Podaj ilość iteracji:");
    std::cin >> iters;
    auto x = solveByJacobi(m, b, iters);
    std::printf("Wynik:\n");
    print(x);
}

void gaussJacobiComparison() {
    auto [a, b] = loadMatrixWithVector("test2.txt");

    int iters = 1;
    puts("Podaj ilość iteracji do metody Jacobiego:");
    std::cin >> iters;

    auto jacobiResult = solveByJacobi(a, b, iters);
    gaussEliminateToUpper(a, b);
    auto gaussResult = solveUpper(a, b);

    std::printf("Wynik Jacobi:\n");
    print(jacobiResult);
    std::printf("Wynik Gauss:\n");
    print(gaussResult);

    for (int i = 0; i < gaussResult.size(); i++) {
        gaussResult[i] -= jacobiResult[i];
        gaussResult[i] = std::fabs(gaussResult[i]);
    }
    std::printf("Błąd bezwględny:\n");
    print(gaussResult);
}

int main() {
    int choose;
    std::puts("Wybierz scenariusz 1-3 (inna wartosc konczy program)");
    std::puts("1 - Gauss");
    std::puts("2 - Jacobi");
    std::puts("3 - Porownanie powyzszych metod");
	std::cin >> choose;
    switch (choose) {
    case 1: gauss();
		break;
	case 2: jacobi();
		break;
	case 3: gaussJacobiComparison();
		break;
	default:
		return 0;
    }
}

