#include "Element.h"
#include "Grid.h"
#include "EquationSystem.h"
#include "GridSelector.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

// Helper function to clean up numerical errors (round very small values to zero)
double cleanValue(double value, double tolerance = 1e-10) {
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

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    // Select grid file
    GridSelector selector("grids");
    std::string gridFile = selector.selectGrid();
    
    if (gridFile.empty()) {
        std::cout << "Exiting program.\n";
        return 0;
    }
    
    // Load the grid
    Grid grid;
    if (!grid.loadFromFile(gridFile)) {
        std::cerr << "Failed to load grid file: " << gridFile << "\n";
        return 1;
    }
    
    std::cout << "\n=== Successfully loaded: " << gridFile << " ===\n";
    
    // Get global data including conductivity
    const auto& globalData = grid.getGlobalData();
    double conductivity = globalData.getConductivity();
    
    const auto& nodes = grid.getNodes();
    const auto& elements = grid.getElements();
    const auto& boundaryConditions = grid.getBoundaryConditions();
    
    // Display boundary conditions
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "=== BOUNDARY CONDITIONS ===" << "\n";
    std::cout << std::string(80, '=') << "\n\n";
    
    std::cout << "Nodes with boundary conditions: ";
    bool first = true;
    for (int nodeId : boundaryConditions) {
        if (!first) std::cout << ", ";
        std::cout << nodeId;
        first = false;
    }
    std::cout << "\n";
    std::cout << "Total boundary nodes: " << boundaryConditions.size() << "\n";
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    std::cout << "Local H and Hbc matrices for each element:\n\n";
    
    // Calculate and display H and Hbc matrices for each element
    for (const auto& element : elements) {
        std::cout << "=== Element " << element.getId() << " ===\n";
        
        // Get node IDs for this element
        const auto& nodeIds = element.getNodeIds();
        std::cout << "Nodes: ";
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            std::cout << nodeIds[i];
            if (i < nodeIds.size() - 1) std::cout << ", ";
        }
        std::cout << "\n\n";
        
        // Extract node coordinates
        std::vector<double> nodeX(4), nodeY(4);
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            const Node& node = nodes[nodeIds[i] - 1];
            nodeX[i] = node.getX();
            nodeY[i] = node.getY();
        }
        
        // Get boundary edges for this element
        auto boundaryEdges = grid.getElementBoundaryEdges(element);
        bool hasBoundary = false;
        for (bool edge : boundaryEdges) {
            if (edge) {
                hasBoundary = true;
                break;
            }
        }
        
        std::cout << "Boundary edges: ";
        if (hasBoundary) {
            std::vector<std::string> edgeNames = {"Bottom", "Right", "Top", "Left"};
            bool first = true;
            for (size_t i = 0; i < boundaryEdges.size(); ++i) {
                if (boundaryEdges[i]) {
                    if (!first) std::cout << ", ";
                    std::cout << edgeNames[i];
                    first = false;
                }
            }
            std::cout << "\n\n";
        } else {
            std::cout << "None\n\n";
        }
        
        try {
            // Calculate and display H matrix
            auto H_matrix = element.calculateHMatrix(nodeX, nodeY, conductivity, 2);
            printMatrix(H_matrix, "Local H Matrix [4x4]");
            
            // Calculate and display Hbc matrix
            double alfa = globalData.getAlfa();
            auto Hbc_matrix = element.calculateHbcMatrix(nodeX, nodeY, alfa, boundaryEdges);
            printMatrix(Hbc_matrix, "Local Hbc Matrix [4x4]");
            
            // Calculate and display P vector
            double tot = globalData.getTot();
            auto P_vector = element.calculatePVector(nodeX, nodeY, alfa, tot, boundaryEdges);
            printVector(P_vector, "Local P Vector [4]");
            
        } catch (const std::exception& e) {
            std::cerr << "Error calculating matrices: " << e.what() << "\n";
        }
        
        std::cout << "\n" << std::string(80, '-') << "\n\n";
    }
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Assemble global equation system
    EquationSystem eqSystem = grid.assembleGlobalEquationSystem(2);
    
    // Print the global H matrix
    eqSystem.printHMatrix();
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Print the global Hbc matrix
    eqSystem.printHbcMatrix();
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Print the global H+Hbc matrix
    eqSystem.printHTotalMatrix();
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Print the global P vector
    eqSystem.printPVector();
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Show element-to-global mapping for verification
    std::cout << "=== ELEMENT TO GLOBAL NODE MAPPING ===\n\n";
    std::cout << "Element | Local Node 1 | Local Node 2 | Local Node 3 | Local Node 4\n";
    std::cout << std::string(70, '-') << "\n";
    
    for (const auto& element : elements) {
        const auto& nodeIds = element.getNodeIds();
        std::cout << std::setw(7) << element.getId() << " | ";
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            std::cout << std::setw(12) << ("Global " + std::to_string(nodeIds[i]));
            if (i < nodeIds.size() - 1) std::cout << " | ";
        }
        std::cout << "\n";
    }
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Verify matrix properties using EquationSystem methods
    eqSystem.printSystemInfo();
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Solve the equation system [H+Hbc]{t} = {P}
    std::cout << "Solving the equation system [H+Hbc]{t} = {P}...\n\n";
    auto temperatures = eqSystem.solve();
    
    // Display the solution
    eqSystem.printSolution(temperatures);
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Export results to CSV for visualization
    std::string outputFile = "results.csv";
    eqSystem.exportResults(outputFile, temperatures, nodes);
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Launch Python visualization
    std::cout << "Launching temperature visualization...\n";
    std::string pythonCmd = "python3 visualize_temperature.py " + outputFile + " " + gridFile;
    int result = system(pythonCmd.c_str());
    
    if (result != 0) {
        std::cerr << "\nNote: Visualization failed. Make sure matplotlib and numpy are installed.\n";
        std::cerr << "You can manually run: python3 visualize_temperature.py " << outputFile << " " << gridFile << "\n";
    }
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    return 0;
}
