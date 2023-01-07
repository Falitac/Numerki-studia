#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

class Matrix {
public:
    Matrix(size_t row, size_t col, double val = 0.f)
        : row(row)
        , col(col)
        , data(row* col, val) {
    }

    double* operator[](size_t n) {
        return &data[n * col];
    }
    size_t rowSize() const { return row; };
    size_t colSize() const { return col; };

    void print() {
        for (int i = 0; i < rowSize(); i++) {
            for (int j = 0; j < colSize(); j++) {
                std::printf("%15.3lf ", (*this)[i][j]);
            }
            std::printf("\n");
        }
    }
private:
    size_t row;
    size_t col;
    std::vector<double> data;
};

void printVec(std::vector<double>& a) {
    for (int i = 0; i < a.size(); i++) {
        std::printf("%15.20lf ", a[i]);
    }
    puts("");
}

int pivotOnDiagonal(Matrix& a, std::vector<double>& b, int diagonal) {
    int highest = diagonal;
    double maxVal = std::fabs(a[diagonal][diagonal]);
    for (int i = diagonal + 1; i < a.rowSize(); i++) {
        auto absVal = std::fabs(a[i][diagonal]);
        if (absVal > maxVal) {
            maxVal = absVal;
            highest = i;
        }
    }
    return highest;
}

void swapRows(Matrix& a, std::vector<double>& b, int rowIndex1, int rowIndex2) {
    if (rowIndex1 == rowIndex2) return;
    std::swap(b[rowIndex1], b[rowIndex2]);
    std::swap_ranges(a[rowIndex1], a[rowIndex1] + a.colSize(), a[rowIndex2]);
}

void gaussEliminateToUpperPivoting(Matrix& a, std::vector<double>& b) {
    int n = a.rowSize();
    for (int i = 0; i < n - 1; i++) {
        int highest = pivotOnDiagonal(a, b, i);
        swapRows(a, b, i, highest);
        for (int j = i + 1; j < n; j++) {
            double m = a[j][i] / a[i][i];
            b[j] -= b[i] * m;
            for (int k = i; k < n; k++) {
                a[j][k] -= a[i][k] * m;
            }
        }
    }
}

std::vector<double> solveUpper(Matrix& a, std::vector<double>& b) {
    int n = a.rowSize();
    std::vector<double> x(n);

    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = n - 1; j > i; j--) {
            x[i] -= a[i][j] * x[j];
        }
        x[i] /= a[i][i];
    }
    return x;
}

struct Point {
    double x;
    double y;
    double w = 1.;
};

std::vector<Point> loadPoints(const std::string& path) {
    std::fstream input(path);
    int n;
    input >> n;

    std::vector<Point> result(n);
    for (int i = 0; i < n; i++) {
        input >> result[i].x >> result[i].y;
        std::printf("%.2f %.2f %.2f\n", result[i].x, result[i].y, result[i].w);
    }

    return result;
}

Matrix getG(std::vector<Point>& points, int n) {
    Matrix g(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < points.size(); k++) {
                g[i][j] += std::pow(points[k].x, (double)i) * std::pow(points[k].x, (double)j) * points[k].w;
            }
        }
    }
    return g;
}

std::vector<double> getF(std::vector<Point>& points, int n) {
    std::vector<double> result(n, 0.);

    for (int i = 0; i < result.size(); i++) {
        for (int j = 0; j < points.size(); j++) {
            result[i] += std::pow(points[j].x, (double)i) * points[j].y * points[j].w;
        }
    }

    return result;
}

int main() {
    auto points = loadPoints("data0.txt");
    int n;
    std::printf("N:");

    std::cin >> n;

    auto g = getG(points, n);
    auto f = getF(points, n);

    puts("G:");
    g.print();

    puts("F:");
    printVec(f);

    gaussEliminateToUpperPivoting(g, f);
    auto a = solveUpper(g, f);

    puts("A:");
    printVec(a);

    return 0;
}