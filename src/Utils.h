#pragma once

#include <vector>
#include <string>

/**
 * @brief Utility functions for numerical computations and output formatting
 */
namespace Utils {

/**
 * @brief Clean up numerical errors (round very small values to zero)
 * @param value The value to clean
 * @param tolerance The tolerance below which values are rounded to zero
 * @return Cleaned value
 */
double cleanValue(double value, double tolerance = 1e-10);

/**
 * @brief Print matrix with formatting
 * @param matrix The matrix to print
 * @param name The name of the matrix
 */
void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name);

/**
 * @brief Print vector with formatting
 * @param vector The vector to print
 * @param name The name of the vector
 */
void printVector(const std::vector<double>& vector, const std::string& name);

} // namespace Utils

#endif // UTILS_H
