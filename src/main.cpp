#include "Element.h"
#include "Grid.h"
#include <iostream>
#include <iomanip>
#include <cmath>

void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name) {
    std::cout << name << ":\n";
    for (const auto& row : matrix) {
        std::cout << "  [";
        for (size_t j = 0; j < row.size(); ++j) {
            std::cout << std::setw(10) << std::fixed << std::setprecision(6) << row[j];
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
    
    std::cout << "=== Shape Function Derivatives and Jacobian Calculations ===\n\n";
    
    // First show the shape function derivative matrices
    std::cout << "2x2 Gauss Points Integration Scheme:\n";
    std::cout << "Gauss points: ±1/sqrt(3) = ±" << 1.0/std::sqrt(3.0) << "\n\n";
    
    // Show Gauss points order and coordinates
    std::cout << "Gauss Points Order (Point Number, ξ, η):\n";
    std::cout << "Point 1: ξ = -1/√3, η = -1/√3  (" << std::setw(8) << -1.0/std::sqrt(3.0) << ", " << std::setw(8) << -1.0/std::sqrt(3.0) << ")\n";
    std::cout << "Point 2: ξ = +1/√3, η = -1/√3  (" << std::setw(8) <<  1.0/std::sqrt(3.0) << ", " << std::setw(8) << -1.0/std::sqrt(3.0) << ")\n";
    std::cout << "Point 3: ξ = -1/√3, η = +1/√3  (" << std::setw(8) << -1.0/std::sqrt(3.0) << ", " << std::setw(8) <<  1.0/std::sqrt(3.0) << ")\n";
    std::cout << "Point 4: ξ = +1/√3, η = +1/√3  (" << std::setw(8) <<  1.0/std::sqrt(3.0) << ", " << std::setw(8) <<  1.0/std::sqrt(3.0) << ")\n\n";
    
    std::vector<std::string> point2Names = {
        "(-1/sqrt(3), -1/sqrt(3))",
        "( 1/sqrt(3), -1/sqrt(3))",
        "(-1/sqrt(3),  1/sqrt(3))",
        "( 1/sqrt(3),  1/sqrt(3))"
    };
    
    // Show dN/dXi matrix
    std::cout << "dN/dXi Matrix (4x4) - Shape function derivatives with respect to Xi:\n";
    auto dN_dXi_2pt = Element::calculateShapeFunctionDerivativeMatrix_Xi(2);
    for (size_t i = 0; i < dN_dXi_2pt.size(); ++i) {
        std::cout << "Point " << (i+1) << " " << point2Names[i] << ": [";
        for (size_t j = 0; j < dN_dXi_2pt[i].size(); ++j) {
            std::cout << std::setw(8) << dN_dXi_2pt[i][j];
            if (j < dN_dXi_2pt[i].size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
    
    // Show dN/dEta matrix
    std::cout << "\ndN/dEta Matrix (4x4) - Shape function derivatives with respect to Eta:\n";
    auto dN_dEta_2pt = Element::calculateShapeFunctionDerivativeMatrix_Eta(2);
    for (size_t i = 0; i < dN_dEta_2pt.size(); ++i) {
        std::cout << "Point " << (i+1) << " " << point2Names[i] << ": [";
        for (size_t j = 0; j < dN_dEta_2pt[i].size(); ++j) {
            std::cout << std::setw(8) << dN_dEta_2pt[i][j];
            if (j < dN_dEta_2pt[i].size() - 1) std::cout << ", ";
        }
        std::cout << " ]\n";
    }
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    const auto& nodes = grid.getNodes();
    const auto& elements = grid.getElements();
    
    std::cout << "Jacobian calculations for elements using 2x2 Gauss points:\n\n";
    
    // Process first 3 elements with 2x2 integration only
    for (int elemIdx = 0; elemIdx < 3 && elemIdx < (int)elements.size(); ++elemIdx) {
        const auto& element = elements[elemIdx];
        
        std::cout << "=== Element " << element.getId() << " ===\n";
        
        // Get node IDs for this element
        const auto& nodeIds = element.getNodeIds();
        std::cout << "Nodes: ";
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            std::cout << nodeIds[i];
            if (i < nodeIds.size() - 1) std::cout << ", ";
        }
        std::cout << "\n";
        
        // Extract node coordinates
        std::vector<double> nodeX(4), nodeY(4);
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            // Find the node with matching ID (nodes are 1-indexed in file, but stored 0-indexed)
            const Node& node = nodes[nodeIds[i] - 1];
            nodeX[i] = node.getX();
            nodeY[i] = node.getY();
        }
        
        std::cout << "Node coordinates:\n";
        for (size_t i = 0; i < 4; ++i) {
            std::cout << "  Node " << nodeIds[i] << ": (" << std::setw(10) << nodeX[i] 
                      << ", " << std::setw(10) << nodeY[i] << ")\n";
        }
        
        // Calculate Jacobians with 2x2 Gauss points
        try {
            auto jacobianData = element.calculateElementJacobians(nodeX, nodeY, 2);
            
            std::cout << "\nJacobian calculations for " << jacobianData.size() << " integration points:\n";
            
            for (size_t i = 0; i < jacobianData.size(); ++i) {
                std::cout << "\nPoint " << (i+1) << " " << point2Names[i] << ":\n";
                printMatrix(jacobianData[i].jacobian, "  Jacobian [J]");
                std::cout << "  Determinant |J| = " << std::setw(12) << jacobianData[i].detJ << "\n";
                printMatrix(jacobianData[i].inverseJacobian, "  Inverse Jacobian [J^-1]");
                
                // Verify that J * J^-1 = I (identity matrix)
                std::cout << "  Verification (J * J^-1):\n";
                for (int row = 0; row < 2; ++row) {
                    std::cout << "  [";
                    for (int col = 0; col < 2; ++col) {
                        double sum = 0.0;
                        for (int k = 0; k < 2; ++k) {
                            sum += jacobianData[i].jacobian[row][k] * jacobianData[i].inverseJacobian[k][col];
                        }
                        std::cout << std::setw(10) << sum;
                        if (col < 1) std::cout << ", ";
                    }
                    std::cout << " ]\n";
                }
            }
            
        } catch (const std::exception& e) {
            std::cerr << "Error calculating Jacobians: " << e.what() << "\n";
        }
        
        std::cout << "\n" << std::string(80, '-') << "\n\n";
    }
    
    // Summary for all elements (2x2 only)
    std::cout << "=== Summary: Jacobian Determinants for All Elements (2x2 Gauss) ===\n\n";
    
    for (const auto& element : elements) {
        const auto& nodeIds = element.getNodeIds();
        
        // Extract node coordinates
        std::vector<double> nodeX(4), nodeY(4);
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            const Node& node = nodes[nodeIds[i] - 1];
            nodeX[i] = node.getX();
            nodeY[i] = node.getY();
        }
        
        try {
            auto jacobianData = element.calculateElementJacobians(nodeX, nodeY, 2);
            
            std::cout << "Element " << std::setw(2) << element.getId() << " Det J = [";
            for (size_t i = 0; i < jacobianData.size(); ++i) {
                std::cout << std::setw(10) << jacobianData[i].detJ;
                if (i < jacobianData.size() - 1) std::cout << ", ";
            }
            std::cout << " ]\n";
            
        } catch (const std::exception& e) {
            std::cerr << "Element " << element.getId() << " - Error: " << e.what() << "\n";
        }
    }
    
    return 0;
}
