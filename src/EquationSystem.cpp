#include "EquationSystem.h"
#include "Node.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <fstream>
#include <omp.h>

// Helper function to clean up numerical errors
static double cleanValue(double value, double tolerance = 1e-10) {
    return (std::abs(value) < tolerance) ? 0.0 : value;
}

EquationSystem::EquationSystem() : systemSize(0) {}

EquationSystem::EquationSystem(int size) : systemSize(size) {
    initializeHMatrix(size);
    initializeHbcMatrix(size);
    initializeCMatrix(size);
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

void EquationSystem::initializeCMatrix(int size) {
    systemSize = size;
    C_global.clear();
    C_global.resize(size, std::vector<double>(size, 0.0));
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

void EquationSystem::setCMatrix(const std::vector<std::vector<double>>& C) {
    C_global = C;
    if (!C_global.empty()) {
        systemSize = C_global.size();
    }
}

void EquationSystem::addToCMatrix(int i, int j, double value) {
    if (i >= 0 && i < systemSize && j >= 0 && j < systemSize) {
        C_global[i][j] += value;
    }
}

void EquationSystem::setPVector(const std::vector<double>& P) {
    P_global = P;
    if (!P_global.empty()) {
        systemSize = P_global.size();
    }
}

void EquationSystem::setPVector(int i, double value) {
    if (i >= 0 && i < systemSize) {
        P_global[i] = value;
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

void EquationSystem::printMatrixHelper(const std::vector<std::vector<double>>& matrix, 
                                       const std::string& name) const {
    if (matrix.empty()) {
        std::cout << name << " empty\n";
        return;
    }
    
    std::cout << name << " [" << matrix.size() << "x" << matrix[0].size() << "]:\n\n";
    
    // Print column headers
    std::cout << "        ";
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        std::cout << std::setw(9) << ("Node " + std::to_string(j+1));
    }
    std::cout << "\n";
    
    // Print matrix with row headers
    for (size_t i = 0; i < matrix.size(); ++i) {
        std::cout << "Node " << std::setw(2) << (i+1) << " ";
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            std::cout << std::setw(9) << std::fixed << std::setprecision(3) << cleanValue(matrix[i][j]);
        }
        std::cout << "\n";
    }
}

void EquationSystem::printHMatrix() const {
    printMatrixHelper(H_global, "H Matrix");
}

void EquationSystem::printHbcMatrix() const {
    printMatrixHelper(Hbc_global, "Hbc Matrix");
}

void EquationSystem::printCMatrix() const {
    printMatrixHelper(C_global, "C Matrix");
}

void EquationSystem::printHTotalMatrix() const {
    auto HTotal = getHTotalMatrix();
    
    if (HTotal.empty()) {
        std::cout << "H+Hbc matrix empty\n";
        return;
    }
    
    std::cout << "H+Hbc Matrix [" << HTotal.size() << "x" << HTotal[0].size() << "]:\n\n";
    
    // Print column headers
    std::cout << "        ";
    for (size_t j = 0; j < HTotal[0].size(); ++j) {
        std::cout << std::setw(9) << ("Node " + std::to_string(j+1));
    }
    std::cout << "\n";
    
    // Print matrix with row headers
    for (size_t i = 0; i < HTotal.size(); ++i) {
        std::cout << "Node " << std::setw(2) << (i+1) << " ";
        for (size_t j = 0; j < HTotal[i].size(); ++j) {
            std::cout << std::setw(9) << std::fixed << std::setprecision(3) << cleanValue(HTotal[i][j]);
        }
        std::cout << "\n";
    }
}

void EquationSystem::printPVector() const {
    if (P_global.empty()) {
        std::cout << "P vector empty\n";
        return;
    }
    
    std::cout << "P Vector [" << P_global.size() << "]:\n\n";
    
    // Print vector with node labels
    for (size_t i = 0; i < P_global.size(); ++i) {
        std::cout << "Node " << std::setw(2) << (i+1) << ": " 
                  << std::setw(12) << std::fixed << std::setprecision(6) << cleanValue(P_global[i]) << "\n";
    }
}

void EquationSystem::printSystemInfo() const {
    std::cout << "SYSTEM INFO\n\n";
    std::cout << "Size: " << systemSize << "x" << systemSize << "\n";
    std::cout << "Symmetric: " << (isSymmetric() ? "Yes" : "No") << "\n";
    std::cout << "Trace: " << std::fixed << std::setprecision(6) << getTrace() << "\n";
    std::cout << "Sum: " << std::fixed << std::setprecision(6) << getSumOfElements() << "\n";
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

std::vector<double> EquationSystem::solveLinearSystem(std::vector<std::vector<double>>& A) const {
    // Solve linear system using Gaussian elimination with partial pivoting
    // Input: augmented matrix [A | b]
    // Output: solution vector x
    
    int n = A.size();
    if (n == 0 || A[0].size() != static_cast<size_t>(n + 1)) {
        std::cerr << "Error: Invalid augmented matrix dimensions.\n";
        return std::vector<double>();
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
            std::cerr << "Warning: Singular matrix at row " << k << "\n";
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
    std::vector<double> solution(n);
    
    for (int i = n - 1; i >= 0; --i) {
        double sum = A[i][n];
        for (int j = i + 1; j < n; ++j) {
            sum -= A[i][j] * solution[j];
        }
        
        if (std::abs(A[i][i]) < 1e-12) {
            std::cerr << "Warning: Zero diagonal at row " << i << "\n";
            solution[i] = 0.0;
        } else {
            solution[i] = sum / A[i][i];
        }
    }
    
    // Clean up numerical errors
    for (int i = 0; i < n; ++i) {
        solution[i] = cleanValue(solution[i]);
    }
    
    return solution;
}

std::vector<double> EquationSystem::solve() const {
    // Solve [H+Hbc]{t} = {P}
    
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
    
    return solveLinearSystem(A);
}

void EquationSystem::printSolution(const std::vector<double>& temperatures) const {
    if (temperatures.empty()) {
        std::cout << "No solution\n";
        return;
    }
    
    std::cout << "TEMPERATURE DISTRIBUTION\n\n";
    std::cout << "Node | Temp [C]\n";
    std::cout << std::string(25, '-') << "\n";
    
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
        
        std::cout << "\n" << std::string(25, '-') << "\n";
        std::cout << "Min: " << std::fixed << std::setprecision(4) 
                  << minTemp << " deg C (node " << (minNode + 1) << ")\n";
        std::cout << "Max: " << std::fixed << std::setprecision(4) 
                  << maxTemp << " deg C (node " << (maxNode + 1) << ")\n";
    }
}

std::vector<std::vector<double>> EquationSystem::solveTransient(double simulationTime, 
                                                                double stepTime,
                                                                double initialTemp) const {
    int n = systemSize;
    int numSteps = static_cast<int>(simulationTime / stepTime);
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "TRANSIENT SIMULATION\n";
    std::cout << std::string(60, '=') << "\n\n";
    std::cout << "Time: " << simulationTime << " s, Step: " << stepTime 
              << " s, Steps: " << numSteps << "\n";
    std::cout << "T0: " << initialTemp << " deg C\n";
    std::cout << "OpenMP threads: " << omp_get_max_threads() << "\n\n";
    
    // Store temperature solutions for each time step
    std::vector<std::vector<double>> temperatureHistory;
    
    // Initialize temperature vector with initial temperature
    std::vector<double> t_current(n, initialTemp);
    temperatureHistory.push_back(t_current);
    
    // Pre-compute [C]/dt + [H+Hbc] (left-hand side matrix)
    std::vector<std::vector<double>> LHS(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            LHS[i][j] = C_global[i][j] / stepTime + H_global[i][j] + Hbc_global[i][j];
        }
    }
    
    // Time stepping loop
    for (int step = 1; step <= numSteps; ++step) {
        double currentTime = step * stepTime;
        
        std::cout << "Step " << step << "/" << numSteps 
                  << " (t=" << currentTime << "s)\n";
        
        // Compute right-hand side: ([C]/dt){t^(i)} + {P}
        std::vector<double> RHS(n, 0.0);
        
        // Parallelize the matrix-vector multiplication
        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                sum += (C_global[i][j] / stepTime) * t_current[j];
            }
            RHS[i] = sum + P_global[i];
        }
        
        // Solve: [LHS]{t^(i+1)} = {RHS}
        // Create augmented matrix [LHS | RHS]
        std::vector<std::vector<double>> A(n, std::vector<double>(n + 1));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                A[i][j] = LHS[i][j];
            }
            A[i][n] = RHS[i];
        }
        
        // Solve using helper method
        std::vector<double> t_next = solveLinearSystem(A);
        
        // Update current temperature and store
        t_current = t_next;
        temperatureHistory.push_back(t_current);
    }
    
    std::cout << "\nComplete.\n";
    std::cout << std::string(60, '=') << "\n\n";
    
    return temperatureHistory;
}

void EquationSystem::exportTransientResults(
    const std::string& filename,
    const std::vector<std::vector<double>>& temperatureHistory,
    const std::vector<Node>& nodes,
    double stepTime) const {
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot create file " << filename << "\n";
        return;
    }
    
    // Write CSV header: node_id, x, y, t_0, t_50, t_100, ...
    file << "node_id,x,y";
    for (size_t step = 0; step < temperatureHistory.size(); ++step) {
        file << ",t_" << std::fixed << std::setprecision(1) << (step * stepTime);
    }
    file << "\n";
    
    // Write data for each node
    for (size_t i = 0; i < nodes.size(); ++i) {
        file << (i + 1) << ","
             << std::fixed << std::setprecision(6) << nodes[i].getX() << ","
             << std::fixed << std::setprecision(6) << nodes[i].getY();
        
        // Write temperature at each time step
        for (const auto& temps : temperatureHistory) {
            if (i < temps.size()) {
                file << "," << std::fixed << std::setprecision(4) << temps[i];
            }
        }
        file << "\n";
    }
    
    file.close();
    std::cout << "Exported transient results to: " << filename << "\n";
}
