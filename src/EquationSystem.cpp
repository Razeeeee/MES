#include "EquationSystem.h"
#include "Node.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

// Helper function to clean up numerical errors
static double cleanValue(double value, double tolerance = 1e-10) {
    return (std::abs(value) < tolerance) ? 0.0 : value;
}

EquationSystem::EquationSystem() : systemSize(0) {}

EquationSystem::EquationSystem(int size) : systemSize(size) {
    initializeHMatrix(size);
    initializeHbcMatrix(size);
    initializePVector(size);
}

void EquationSystem::initializeHMatrix(int size) {
    systemSize = size;
    H_global.clear();
    H_global.resize(size, std::vector<double>(size, 0.0));
}

void EquationSystem::initializeHbcMatrix(int size) {
    systemSize = size;
    Hbc_global.clear();
    Hbc_global.resize(size, std::vector<double>(size, 0.0));
}

void EquationSystem::initializePVector(int size) {
    systemSize = size;
    P_global.clear();
    P_global.resize(size, 0.0);
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

void EquationSystem::setHbcMatrix(const std::vector<std::vector<double>>& Hbc) {
    Hbc_global = Hbc;
    if (!Hbc_global.empty()) {
        systemSize = Hbc_global.size();
    }
}

void EquationSystem::addToHbcMatrix(int i, int j, double value) {
    if (i >= 0 && i < systemSize && j >= 0 && j < systemSize) {
        Hbc_global[i][j] += value;
    }
}

void EquationSystem::setPVector(const std::vector<double>& P) {
    P_global = P;
    if (!P_global.empty()) {
        systemSize = P_global.size();
    }
}

void EquationSystem::addToPVector(int i, double value) {
    if (i >= 0 && i < systemSize) {
        P_global[i] += value;
    }
}

std::vector<std::vector<double>> EquationSystem::getHTotalMatrix() const {
    std::vector<std::vector<double>> HTotal(systemSize, std::vector<double>(systemSize, 0.0));
    
    for (int i = 0; i < systemSize; ++i) {
        for (int j = 0; j < systemSize; ++j) {
            HTotal[i][j] = H_global[i][j] + Hbc_global[i][j];
        }
    }
    
    return HTotal;
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
            std::cout << std::setw(8) << std::fixed << std::setprecision(3) << cleanValue(H_global[i][j]);
        }
        std::cout << "\n";
    }
}

void EquationSystem::printHbcMatrix() const {
    if (Hbc_global.empty()) {
        std::cout << "Hbc matrix is empty.\n";
        return;
    }
    
    std::cout << "Global Hbc Matrix [" << Hbc_global.size() << " x " << Hbc_global[0].size() << "]:\n\n";
    
    // Print column headers
    std::cout << "        ";
    for (size_t j = 0; j < Hbc_global[0].size(); ++j) {
        std::cout << std::setw(8) << ("Node " + std::to_string(j+1));
    }
    std::cout << "\n";
    
    // Print matrix with row headers
    for (size_t i = 0; i < Hbc_global.size(); ++i) {
        std::cout << "Node " << std::setw(2) << (i+1) << " ";
        for (size_t j = 0; j < Hbc_global[i].size(); ++j) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(3) << cleanValue(Hbc_global[i][j]);
        }
        std::cout << "\n";
    }
}

void EquationSystem::printHTotalMatrix() const {
    auto HTotal = getHTotalMatrix();
    
    if (HTotal.empty()) {
        std::cout << "H+Hbc matrix is empty.\n";
        return;
    }
    
    std::cout << "Global H+Hbc Matrix [" << HTotal.size() << " x " << HTotal[0].size() << "]:\n\n";
    
    // Print column headers
    std::cout << "        ";
    for (size_t j = 0; j < HTotal[0].size(); ++j) {
        std::cout << std::setw(8) << ("Node " + std::to_string(j+1));
    }
    std::cout << "\n";
    
    // Print matrix with row headers
    for (size_t i = 0; i < HTotal.size(); ++i) {
        std::cout << "Node " << std::setw(2) << (i+1) << " ";
        for (size_t j = 0; j < HTotal[i].size(); ++j) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(3) << cleanValue(HTotal[i][j]);
        }
        std::cout << "\n";
    }
}

void EquationSystem::printPVector() const {
    if (P_global.empty()) {
        std::cout << "P vector is empty.\n";
        return;
    }
    
    std::cout << "Global P Vector [" << P_global.size() << "]:\n\n";
    
    // Print vector with node labels
    for (size_t i = 0; i < P_global.size(); ++i) {
        std::cout << "Node " << std::setw(2) << (i+1) << ": " 
                  << std::setw(12) << std::fixed << std::setprecision(6) << cleanValue(P_global[i]) << "\n";
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

std::vector<double> EquationSystem::solve() const {
    // Solve [H+Hbc]{t} = {P}
    // Using Gaussian elimination with partial pivoting
    
    int n = systemSize;
    
    // Create augmented matrix [H+Hbc | P]
    std::vector<std::vector<double>> A(n, std::vector<double>(n + 1));
    
    // Fill the augmented matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = H_global[i][j] + Hbc_global[i][j];
        }
        A[i][n] = P_global[i]; // Right-hand side (P vector)
    }
    
    // Forward elimination with partial pivoting
    for (int k = 0; k < n - 1; ++k) {
        // Find pivot
        int pivotRow = k;
        double maxVal = std::abs(A[k][k]);
        
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(A[i][k]) > maxVal) {
                maxVal = std::abs(A[i][k]);
                pivotRow = i;
            }
        }
        
        // Swap rows if necessary
        if (pivotRow != k) {
            std::swap(A[k], A[pivotRow]);
        }
        
        // Check for singular matrix
        if (std::abs(A[k][k]) < 1e-12) {
            std::cerr << "Warning: Matrix is singular or nearly singular at row " << k << "\n";
            continue;
        }
        
        // Eliminate column k
        for (int i = k + 1; i < n; ++i) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j <= n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
        }
    }
    
    // Back substitution
    std::vector<double> temperatures(n);
    
    for (int i = n - 1; i >= 0; --i) {
        double sum = A[i][n];
        for (int j = i + 1; j < n; ++j) {
            sum -= A[i][j] * temperatures[j];
        }
        
        if (std::abs(A[i][i]) < 1e-12) {
            std::cerr << "Warning: Zero diagonal element at row " << i << "\n";
            temperatures[i] = 0.0;
        } else {
            temperatures[i] = sum / A[i][i];
        }
    }
    
    // Clean up numerical errors
    for (int i = 0; i < n; ++i) {
        temperatures[i] = cleanValue(temperatures[i]);
    }
    
    return temperatures;
}

void EquationSystem::printSolution(const std::vector<double>& temperatures) const {
    if (temperatures.empty()) {
        std::cout << "Solution vector is empty.\n";
        return;
    }
    
    std::cout << "=== TEMPERATURE DISTRIBUTION ===\n\n";
    std::cout << "Node |  Temperature [°C]\n";
    std::cout << std::string(30, '-') << "\n";
    
    for (size_t i = 0; i < temperatures.size(); ++i) {
        std::cout << std::setw(4) << (i + 1) << " | " 
                  << std::setw(16) << std::fixed << std::setprecision(6) 
                  << temperatures[i] << "\n";
    }
    
    // Find min and max temperatures
    if (!temperatures.empty()) {
        double minTemp = temperatures[0];
        double maxTemp = temperatures[0];
        int minNode = 0, maxNode = 0;
        
        for (size_t i = 1; i < temperatures.size(); ++i) {
            if (temperatures[i] < minTemp) {
                minTemp = temperatures[i];
                minNode = i;
            }
            if (temperatures[i] > maxTemp) {
                maxTemp = temperatures[i];
                maxNode = i;
            }
        }
        
        std::cout << "\n" << std::string(30, '-') << "\n";
        std::cout << "Minimum temperature: " << std::fixed << std::setprecision(6) 
                  << minTemp << " °C at node " << (minNode + 1) << "\n";
        std::cout << "Maximum temperature: " << std::fixed << std::setprecision(6) 
                  << maxTemp << " °C at node " << (maxNode + 1) << "\n";
    }
}

void EquationSystem::exportResults(const std::string& filename, 
                                   const std::vector<double>& temperatures,
                                   const std::vector<Node>& nodes) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }
    
    // Write CSV header
    file << "node_id,x,y,temperature\n";
    
    // Write data for each node
    for (size_t i = 0; i < temperatures.size() && i < nodes.size(); ++i) {
        file << (i + 1) << ","
             << std::fixed << std::setprecision(10) << nodes[i].getX() << ","
             << std::fixed << std::setprecision(10) << nodes[i].getY() << ","
             << std::fixed << std::setprecision(6) << temperatures[i] << "\n";
    }
    
    file.close();
    std::cout << "Results exported to: " << filename << "\n";
}
