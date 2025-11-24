#pragma once
#include <vector>
#include <string>

// Forward declaration
class Node;

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
    std::vector<std::vector<double>> Hbc_global; // Global Hbc matrix [N x N]
    std::vector<double> P_global;                // Global P vector [N]
    int systemSize;                              // Number of nodes (N)

public:
    EquationSystem();
    explicit EquationSystem(int size);
    
    // Matrix operations
    void setHMatrix(const std::vector<std::vector<double>>& H);
    void addToHMatrix(int i, int j, double value);
    void initializeHMatrix(int size);
    
    void setHbcMatrix(const std::vector<std::vector<double>>& Hbc);
    void addToHbcMatrix(int i, int j, double value);
    void initializeHbcMatrix(int size);
    
    // P vector operations
    void setPVector(const std::vector<double>& P);
    void addToPVector(int i, double value);
    void initializePVector(int size);
    
    // Getters
    const std::vector<std::vector<double>>& getHMatrix() const { return H_global; }
    const std::vector<std::vector<double>>& getHbcMatrix() const { return Hbc_global; }
    const std::vector<double>& getPVector() const { return P_global; }
    std::vector<std::vector<double>> getHTotalMatrix() const; // Returns H + Hbc
    int getSystemSize() const { return systemSize; }
    
    // Display methods
    void printHMatrix() const;
    void printHbcMatrix() const;
    void printHTotalMatrix() const; // Print H + Hbc
    void printPVector() const; // Print P vector
    void printSystemInfo() const;
    
    // Solver methods
    // Solve the equation system [H+Hbc]{t} = -{P}
    // Returns temperature vector {t}
    std::vector<double> solve() const;
    void printSolution(const std::vector<double>& temperatures) const;
    
    // Matrix properties
    bool isSymmetric(double tolerance = 1e-10) const;
    double getTrace() const;
    double getSumOfElements() const;
    
    // Export methods
    void exportResults(const std::string& filename, 
                      const std::vector<double>& temperatures,
                      const std::vector<Node>& nodes) const;
};
