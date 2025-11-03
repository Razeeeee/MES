#include "EquationSystem.h"
#include <iostream>
#include <iomanip>
#include <cmath>

EquationSystem::EquationSystem() : systemSize(0) {}

EquationSystem::EquationSystem(int size) : systemSize(size) {
    initializeHMatrix(size);
}

void EquationSystem::initializeHMatrix(int size) {
    systemSize = size;
    H_global.clear();
    H_global.resize(size, std::vector<double>(size, 0.0));
}

void EquationSystem::setHMatrix(const std::vector<std::vector<double>>& H) {
    H_global = H;
    if (!H_global.empty()) {
        systemSize = H_global.size();
    }
}

void EquationSystem::addToHMatrix(int i, int j, double value) {
    if (i >= 0 && i < systemSize && j >= 0 && j < systemSize) {
        H_global[i][j] += value;
    }
}

void EquationSystem::printHMatrix() const {
    if (H_global.empty()) {
        std::cout << "H matrix is empty.\n";
        return;
    }
    
    std::cout << "Global H Matrix [" << H_global.size() << " x " << H_global[0].size() << "]:\n\n";
    
    // Print column headers
    std::cout << "        ";
    for (size_t j = 0; j < H_global[0].size(); ++j) {
        std::cout << std::setw(8) << ("Node " + std::to_string(j+1));
    }
    std::cout << "\n";
    
    // Print matrix with row headers
    for (size_t i = 0; i < H_global.size(); ++i) {
        std::cout << "Node " << std::setw(2) << (i+1) << " ";
        for (size_t j = 0; j < H_global[i].size(); ++j) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(3) << H_global[i][j];
        }
        std::cout << "\n";
    }
}

void EquationSystem::printSystemInfo() const {
    std::cout << "=== EQUATION SYSTEM PROPERTIES ===\n\n";
    std::cout << "System Size: " << systemSize << " x " << systemSize << "\n";
    std::cout << "Is Symmetric: " << (isSymmetric() ? "Yes" : "No") << "\n";
    std::cout << "Trace (sum of diagonal): " << std::fixed << std::setprecision(6) << getTrace() << "\n";
    std::cout << "Sum of all elements: " << std::fixed << std::setprecision(6) << getSumOfElements() << "\n";
}

bool EquationSystem::isSymmetric(double tolerance) const {
    if (H_global.empty()) return true;
    
    for (size_t i = 0; i < H_global.size(); ++i) {
        for (size_t j = i + 1; j < H_global[i].size(); ++j) {
            if (std::abs(H_global[i][j] - H_global[j][i]) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

double EquationSystem::getTrace() const {
    double trace = 0.0;
    for (size_t i = 0; i < H_global.size(); ++i) {
        trace += H_global[i][i];
    }
    return trace;
}

double EquationSystem::getSumOfElements() const {
    double sum = 0.0;
    for (const auto& row : H_global) {
        for (double val : row) {
            sum += val;
        }
    }
    return sum;
}
