#pragma once
#include <vector>
#include <string>
#include <set>

// Global FEM equation system with solver methods
class EquationSystem {
private:
    std::vector<std::vector<double>> H_global;
    std::vector<std::vector<double>> Hbc_global;
    std::vector<std::vector<double>> C_global;
    std::vector<double> P_global;
    int systemSize;

    void printMatrixHelper(const std::vector<std::vector<double>>& matrix, 
                          const std::string& name) const;
    std::vector<double> solveLinearSystem(std::vector<std::vector<double>>& augmentedMatrix) const;

public:
    EquationSystem();
    explicit EquationSystem(int size);
    
    void setHMatrix(const std::vector<std::vector<double>>& H);
    void addToHMatrix(int i, int j, double value);
    void initializeHMatrix(int size);
    
    void setHbcMatrix(const std::vector<std::vector<double>>& Hbc);
    void addToHbcMatrix(int i, int j, double value);
    void initializeHbcMatrix(int size);
    
    void setCMatrix(const std::vector<std::vector<double>>& C);
    void addToCMatrix(int i, int j, double value);
    void initializeCMatrix(int size);
    
    void setPVector(const std::vector<double>& P);
    void setPVector(int i, double value);  // Set individual element
    void addToPVector(int i, double value);
    void initializePVector(int size);
    
    const std::vector<std::vector<double>>& getHMatrix() const noexcept { return H_global; }
    const std::vector<std::vector<double>>& getHbcMatrix() const noexcept { return Hbc_global; }
    const std::vector<std::vector<double>>& getCMatrix() const noexcept { return C_global; }
    const std::vector<double>& getPVector() const noexcept { return P_global; }
    std::vector<std::vector<double>> getHTotalMatrix() const;
    int getSystemSize() const noexcept { return systemSize; }
    
    void printHMatrix() const;
    void printHbcMatrix() const;
    void printCMatrix() const;
    void printHTotalMatrix() const;
    void printPVector() const;
    void printSystemInfo() const;
    
    std::vector<double> solve() const;
    std::vector<std::vector<double>> solveTransient(double simulationTime, 
                                                    double stepTime,
                                                    double initialTemp) const;
    void printSolution(const std::vector<double>& temperatures) const;
    
    // Export results to CSV
    void exportTransientResults(const std::string& filename,
                               const std::vector<std::vector<double>>& temperatureHistory,
                               const std::vector<class Node>& nodes,
                               double stepTime) const;
    
    // Calculate power needed to maintain temperature (heat loss through boundary)
    double calculatePower(const std::vector<double>& temperatures,
                         const std::vector<class Node>& nodes,
                         const std::vector<class Element>& elements,
                         const std::set<int>& boundaryConditions,
                         double alfa,
                         double tot) const;
    
    bool isSymmetric(double tolerance = 1e-10) const;
    double getTrace() const;
    double getSumOfElements() const;
};
