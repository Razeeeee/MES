#include "Element.h"
#include "Grid.h"
#include "EquationSystem.h"
#include <iostream>
#include <iomanip>
#include <cmath>

void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name) {
    std::cout << name << ":\n";
    for (const auto& row : matrix) {
        std::cout << "  [";
        for (size_t j = 0; j < row.size(); ++j) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(6) << row[j];
            if (j < row.size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    // Load the grid
    Grid grid;
    if (!grid.loadFromFile("grids/Test2_4_4_MixGrid.txt")) {
        std::cerr << "Failed to load grid file!\n";
        return 1;
    }
    
    // Get global data including conductivity
    const auto& globalData = grid.getGlobalData();
    double conductivity = globalData.getConductivity();
    
    const auto& nodes = grid.getNodes();
    const auto& elements = grid.getElements();
    
    std::cout << "Local H matrices for each element:\n\n";
    
    // Calculate and display H matrix for each element
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
        
        try {
            // Calculate and display H matrix
            auto H_matrix = element.calculateHMatrix(nodeX, nodeY, conductivity, 2);
            printMatrix(H_matrix, "Local H Matrix [4x4]");
            
        } catch (const std::exception& e) {
            std::cerr << "Error calculating H matrix: " << e.what() << "\n";
        }
        
        std::cout << "\n" << std::string(80, '-') << "\n\n";
    }
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Assemble global equation system
    EquationSystem eqSystem = grid.assembleGlobalEquationSystem(2);
    
    // Print the global H matrix
    eqSystem.printHMatrix();
    
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
    
    return 0;
}
