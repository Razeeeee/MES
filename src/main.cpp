#include "Element.h"
#include "Grid.h"
#include "EquationSystem.h"
#include "GridSelector.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <limits>

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
            
            // Calculate and display C matrix
            double density = globalData.getDensity();
            double specificHeat = globalData.getSpecificHeat();
            auto C_matrix = element.calculateCMatrix(nodeX, nodeY, density, specificHeat, 2);
            printMatrix(C_matrix, "Local C Matrix [4x4]");
            
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
    
    // Print the global C matrix
    eqSystem.printCMatrix();
    
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
    
    // =========================================================
    // TRANSIENT ANALYSIS (Time-dependent solution)
    // =========================================================
    
    std::cout << "=== TRANSIENT HEAT TRANSFER ANALYSIS ===\n\n";
    
    // Get time parameters from global data
    double simulationTime = globalData.getSimulationTime();
    double stepTime = globalData.getSimulationStepTime();
    double initialTemp = globalData.getInitialTemp();
    
    // Solve transient heat transfer equation
    auto temperatureHistory = eqSystem.solveTransient(simulationTime, stepTime, initialTemp);
    
    // Display final temperature distribution
    std::cout << "\n=== FINAL TEMPERATURE DISTRIBUTION (t = " << simulationTime << " s) ===\n\n";
    eqSystem.printSolution(temperatureHistory.back());
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Display temperature at intermediate time steps
    std::cout << "=== TEMPERATURE HISTORY ===\n\n";
    int numSteps = temperatureHistory.size();
    std::cout << "Total time steps: " << numSteps << "\n\n";
    
    // Track min and max temperatures throughout the entire simulation
    double minTemp = std::numeric_limits<double>::max();
    double maxTemp = std::numeric_limits<double>::lowest();
    int minTempNode = -1;
    int maxTempNode = -1;
    double minTempTime = 0.0;
    double maxTempTime = 0.0;
    int minTempElement = -1;
    int maxTempElement = -1;
    
    // Analyze all time steps to find global min/max
    for (int step = 0; step < numSteps; ++step) {
        double time = step * stepTime;
        const auto& temps = temperatureHistory[step];
        
        for (size_t i = 0; i < temps.size(); ++i) {
            if (temps[i] < minTemp) {
                minTemp = temps[i];
                minTempNode = i + 1;  // Node IDs are 1-based
                minTempTime = time;
            }
            if (temps[i] > maxTemp) {
                maxTemp = temps[i];
                maxTempNode = i + 1;  // Node IDs are 1-based
                maxTempTime = time;
            }
        }
    }
    
    // Find which element contains the min/max temperature nodes
    for (const auto& element : elements) {
        const auto& nodeIds = element.getNodeIds();
        for (int nodeId : nodeIds) {
            if (nodeId == minTempNode) {
                minTempElement = element.getId();
            }
            if (nodeId == maxTempNode) {
                maxTempElement = element.getId();
            }
        }
    }
    
    // Display global min/max temperature information
    std::cout << "=== GLOBAL TEMPERATURE EXTREMES (ENTIRE SIMULATION) ===\n\n";
    std::cout << "Minimum Temperature:\n";
    std::cout << "  Temperature: " << std::fixed << std::setprecision(4) << minTemp << " °C\n";
    std::cout << "  Node:        " << minTempNode << "\n";
    std::cout << "  Element:     " << minTempElement << "\n";
    std::cout << "  Time:        " << std::fixed << std::setprecision(1) << minTempTime << " s\n\n";
    
    std::cout << "Maximum Temperature:\n";
    std::cout << "  Temperature: " << std::fixed << std::setprecision(4) << maxTemp << " °C\n";
    std::cout << "  Node:        " << maxTempNode << "\n";
    std::cout << "  Element:     " << maxTempElement << "\n";
    std::cout << "  Time:        " << std::fixed << std::setprecision(1) << maxTempTime << " s\n\n";
    
    std::cout << std::string(80, '-') << "\n\n";
    
    // Show temperatures at selected nodes for each time step
    std::vector<int> displayNodes = {1, 4, 13, 16}; // Corner nodes
    std::cout << "Time [s] | ";
    for (int nodeId : displayNodes) {
        std::cout << "Node " << std::setw(2) << nodeId << " [°C] | ";
    }
    std::cout << "\n" << std::string(80, '-') << "\n";
    
    for (int step = 0; step < numSteps; ++step) {
        double time = step * stepTime;
        std::cout << std::setw(8) << std::fixed << std::setprecision(1) << time << " | ";
        for (int nodeId : displayNodes) {
            if (nodeId - 1 < static_cast<int>(temperatureHistory[step].size())) {
                std::cout << std::setw(13) << std::fixed << std::setprecision(4) 
                         << temperatureHistory[step][nodeId - 1] << " | ";
            }
        }
        std::cout << "\n";
    }
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // =========================================================
    // STEADY-STATE ANALYSIS (for comparison)
    // =========================================================
    
    // Solve the equation system [H+Hbc]{t} = {P}
    std::cout << "Solving the steady-state equation system [H+Hbc]{t} = {P}...\n\n";
    auto temperatures = eqSystem.solve();
    
    // Display the solution
    std::cout << "=== STEADY-STATE TEMPERATURE DISTRIBUTION ===\n\n";
    eqSystem.printSolution(temperatures);
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Export results to CSV for visualization
    std::string outputFile = "results.csv";
    eqSystem.exportResults(outputFile, temperatures, nodes);
    
    // Export transient results to CSV
    std::string transientOutputFile = "results_transient.csv";
    std::ofstream transientFile(transientOutputFile);
    if (transientFile.is_open()) {
        // Write header with time steps
        transientFile << "node_id,x,y";
        for (size_t step = 0; step < temperatureHistory.size(); ++step) {
            transientFile << ",t_" << (step * stepTime);
        }
        transientFile << "\n";
        
        // Write data for each node
        for (size_t i = 0; i < nodes.size(); ++i) {
            transientFile << (i + 1) << ","
                         << std::fixed << std::setprecision(10) << nodes[i].getX() << ","
                         << std::fixed << std::setprecision(10) << nodes[i].getY();
            
            for (const auto& temps : temperatureHistory) {
                if (i < temps.size()) {
                    transientFile << "," << std::fixed << std::setprecision(6) << temps[i];
                }
            }
            transientFile << "\n";
        }
        transientFile.close();
        std::cout << "Transient results exported to: " << transientOutputFile << "\n";
    }
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    // Launch Python visualization for steady-state
    std::cout << "Launching steady-state temperature visualization...\n";
    std::string pythonCmd = "python3 visualize_temperature.py " + outputFile + " " + gridFile;
    int result = system(pythonCmd.c_str());
    
    if (result != 0) {
        std::cerr << "\nNote: Steady-state visualization failed. Make sure matplotlib and numpy are installed.\n";
        std::cerr << "You can manually run: python3 visualize_temperature.py " << outputFile << " " << gridFile << "\n";
    }
    
    // Launch Python visualization for transient analysis
    std::cout << "\nLaunching transient temperature visualization with time slider...\n";
    std::string pythonCmdTransient = "python3 visualize_temperature_transient.py " + transientOutputFile + " " + gridFile;
    int resultTransient = system(pythonCmdTransient.c_str());
    
    if (resultTransient != 0) {
        std::cerr << "\nNote: Transient visualization failed. Make sure matplotlib and numpy are installed.\n";
        std::cerr << "You can manually run: python3 visualize_temperature_transient.py " << transientOutputFile << " " << gridFile << "\n";
    }
    
    std::cout << "\n" << std::string(80, '=') << "\n\n";
    
    return 0;
}
