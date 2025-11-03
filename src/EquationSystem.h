#pragma once
#include <vector>
#include <string>

/**
 * @brief Class representing a system of equations
 * 
 * This class manages the global matrices and vectors including:
 * - Global H matrix (stiffness/conductivity matrix)
 * - Global P vector (load vector) - to be added
 * - Solution methods - to be added
 */
class EquationSystem {
private:
    std::vector<std::vector<double>> H_global;  // Global H matrix [N x N]
    // std::vector<double> P_global;             // Global P vector [N] - to be added
    int systemSize;                              // Number of nodes (N)

public:
    EquationSystem();
    explicit EquationSystem(int size);
    
    // Matrix operations
    void setHMatrix(const std::vector<std::vector<double>>& H);
    void addToHMatrix(int i, int j, double value);
    void initializeHMatrix(int size);
    
    // Getters
    const std::vector<std::vector<double>>& getHMatrix() const { return H_global; }
    int getSystemSize() const { return systemSize; }
    
    // Display methods
    void printHMatrix() const;
    void printSystemInfo() const;
    
    // Matrix properties
    bool isSymmetric(double tolerance = 1e-10) const;
    double getTrace() const;
    double getSumOfElements() const;
};
