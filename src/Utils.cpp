#include "Utils.h"
#include <iostream>
#include <iomanip>
#include <cmath>

namespace Utils {

double cleanValue(double value, double tolerance) {
    return (std::abs(value) < tolerance) ? 0.0 : value;
}

void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name) {
    std::cout << name << ":\n";
    for (const auto& row : matrix) {
        std::cout << "  [";
        for (size_t j = 0; j < row.size(); ++j) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(6) << cleanValue(row[j]);
            if (j < row.size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
}

void printVector(const std::vector<double>& vector, const std::string& name) {
    std::cout << name << ":\n  [";
    for (size_t i = 0; i < vector.size(); ++i) {
        std::cout << std::setw(12) << std::fixed << std::setprecision(6) << cleanValue(vector[i]);
        if (i < vector.size() - 1) std::cout << ", ";
    }
    std::cout << " ]\n";
}

} // namespace Utils
