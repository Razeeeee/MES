#include "Element.h"
#include "Grid.h"
#include "EquationSystem.h"
#include "GridSelector.h"
#include "Config.h"
#include "Utils.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    Config config("config/fem_config.txt");
    config.print();
    
    GridSelector selector("grids");
    std::string gridFile = selector.selectGrid();
    
    if (gridFile.empty()) {
        std::cout << "Exiting program.\n";
        return 0;
    }
    
    Grid grid;
    if (!grid.loadFromFile(gridFile)) {
        std::cerr << "Error: Failed to load " << gridFile << "\n";
        return 1;
    }
    
    std::cout << "\nLoaded: " << gridFile << "\n";
    
    const auto& globalData = grid.getGlobalData();
    const auto& nodes = grid.getNodes();
    const auto& elements = grid.getElements();
    const auto& boundaryConditions = grid.getBoundaryConditions();
    
    double conductivity = globalData.getConductivity();
    double alfa = globalData.getAlfa();
    double tot = globalData.getTot();
    double density = globalData.getDensity();
    double specificHeat = globalData.getSpecificHeat();
    
    if (config.shouldPrintBoundaryConditions()) {
        std::cout << "\n" << std::string(60, '=') << "\n";
        std::cout << "BOUNDARY CONDITIONS\n";
        std::cout << std::string(60, '=') << "\n\n";
        
        std::cout << "Convection BC nodes: ";
        bool first = true;
        for (int nodeId : boundaryConditions) {
            if (!first) std::cout << ", ";
            std::cout << nodeId;
            first = false;
        }
        std::cout << "\nTotal: " << boundaryConditions.size() << " nodes\n";
        std::cout << std::string(60, '=') << "\n\n";
    }
    
    EquationSystem eqSystem = grid.assembleGlobalEquationSystem(
        config.getGaussPointsH(),
        config.getGaussPointsC(),
        config.getGaussPointsBoundary());
    
    if (config.shouldPrintLocalMatrices()) {
        std::cout << "LOCAL ELEMENT MATRICES\n\n";
        
        for (const auto& element : elements) {
            std::cout << "Element " << element.getId() << " (Nodes: ";
            const auto& nodeIds = element.getNodeIds();
            for (size_t i = 0; i < nodeIds.size(); ++i) {
                std::cout << nodeIds[i];
                if (i < nodeIds.size() - 1) std::cout << ", ";
            }
            std::cout << ")\n";
            
            auto boundaryEdges = grid.getElementBoundaryEdges(element);
            bool hasBoundary = false;
            for (bool edge : boundaryEdges) {
                if (edge) { hasBoundary = true; break; }
            }
            
            if (hasBoundary) {
                std::cout << "Boundary edges: ";
                std::vector<std::string> edgeNames = {"Bottom", "Right", "Top", "Left"};
                bool first = true;
                for (size_t i = 0; i < boundaryEdges.size(); ++i) {
                    if (boundaryEdges[i]) {
                        if (!first) std::cout << ", ";
                        std::cout << edgeNames[i];
                        first = false;
                    }
                }
                std::cout << "\n";
            } else {
                std::cout << "Boundary edges: None\n";
            }
            std::cout << "\n";
            
            Utils::printMatrix(element.getHLocal(), "Local H Matrix [4x4]");
            Utils::printMatrix(element.getCLocal(), "Local C Matrix [4x4]");
            Utils::printMatrix(element.getHbcLocal(), "Local Hbc Matrix [4x4]");
            Utils::printVector(element.getPLocal(), "Local P Vector [4]");
            
            std::cout << "\n" << std::string(60, '-') << "\n\n";
        }
    }
    
    if (config.shouldPrintGlobalMatrices()) {
        std::cout << "\n" << std::string(60, '=') << "\n";
        eqSystem.printHMatrix();
        
        std::cout << "\n" << std::string(60, '=') << "\n\n";
        eqSystem.printHbcMatrix();
        
        std::cout << "\n" << std::string(60, '=') << "\n\n";
        eqSystem.printCMatrix();
        
        std::cout << "\n" << std::string(60, '=') << "\n\n";
        eqSystem.printHTotalMatrix();
        
        std::cout << "\n" << std::string(60, '=') << "\n\n";
        eqSystem.printPVector();
        
        std::cout << "\n" << std::string(60, '=') << "\n\n";
    }
    
    if (config.shouldPrintAssemblyMapping()) {
        std::cout << "ELEMENT-TO-NODE MAPPING\n\n";
        std::cout << "Element | Node 1 | Node 2 | Node 3 | Node 4\n";
        std::cout << std::string(50, '-') << "\n";
        
        for (const auto& element : elements) {
            const auto& nodeIds = element.getNodeIds();
            std::cout << std::setw(7) << element.getId() << " | ";
            for (size_t i = 0; i < nodeIds.size(); ++i) {
                std::cout << std::setw(6) << nodeIds[i];
                if (i < nodeIds.size() - 1) std::cout << " | ";
            }
            std::cout << "\n";
        }
        std::cout << std::string(50, '=') << "\n\n";
    }
    
    if (config.shouldPrintSystemInfo()) {
        eqSystem.printSystemInfo();
        std::cout << "\n" << std::string(60, '=') << "\n\n";
    }
    
    std::cout << "TRANSIENT ANALYSIS\n\n";
    
    double simulationTime = globalData.getSimulationTime();
    double stepTime = globalData.getSimulationStepTime();
    double initialTemp = globalData.getInitialTemp();
    
    auto temperatureHistory = eqSystem.solveTransient(simulationTime, stepTime, initialTemp);
    
    if (config.shouldPrintFinalResults()) {
        std::cout << "\nFINAL RESULTS (t = " << simulationTime << " s)\n\n";
        eqSystem.printSolution(temperatureHistory.back());
        std::cout << "\n" << std::string(60, '=') << "\n\n";
    }
    
    if (config.shouldPrintTimeSteps()) {
        std::cout << "TEMPERATURE HISTORY\n\n";
        int numSteps = temperatureHistory.size();
        std::cout << "Time steps: " << numSteps << "\n\n";
        
        double minTemp = std::numeric_limits<double>::max();
        double maxTemp = std::numeric_limits<double>::lowest();
        int minTempNode = -1, maxTempNode = -1;
        double minTempTime = 0.0, maxTempTime = 0.0;
        int minTempElement = -1, maxTempElement = -1;
        
        for (int step = 0; step < numSteps; ++step) {
            double time = step * stepTime;
            const auto& temps = temperatureHistory[step];
            
            for (size_t i = 0; i < temps.size(); ++i) {
                if (temps[i] < minTemp) {
                    minTemp = temps[i];
                    minTempNode = i + 1;
                    minTempTime = time;
                }
                if (temps[i] > maxTemp) {
                    maxTemp = temps[i];
                    maxTempNode = i + 1;
                    maxTempTime = time;
                }
            }
        }
        
        for (const auto& element : elements) {
            const auto& nodeIds = element.getNodeIds();
            for (int nodeId : nodeIds) {
                if (nodeId == minTempNode) minTempElement = element.getId();
                if (nodeId == maxTempNode) maxTempElement = element.getId();
            }
        }
        
        std::cout << "TEMPERATURE EXTREMES\n\n";
        std::cout << "Min: " << std::fixed << std::setprecision(4) << minTemp << " deg C "
                  << "(Node " << minTempNode << ", Element " << minTempElement 
                  << ", t=" << std::setprecision(1) << minTempTime << "s)\n";
        
        std::cout << "Max: " << std::fixed << std::setprecision(4) << maxTemp << " deg C "
                  << "(Node " << maxTempNode << ", Element " << maxTempElement 
                  << ", t=" << std::setprecision(1) << maxTempTime << "s)\n\n";
        
        std::cout << std::string(60, '-') << "\n\n";
        
        std::cout << "Time [s] | Min Temp [C] (Node) | Max Temp [C] (Node) |\n";
        std::cout << std::string(60, '-') << "\n";
        
        for (int step = 0; step < numSteps; ++step) {
            double time = step * stepTime;
            const auto& temps = temperatureHistory[step];
            
            double stepMin = std::numeric_limits<double>::max();
            double stepMax = std::numeric_limits<double>::lowest();
            int stepMinNode = -1, stepMaxNode = -1;
            
            for (size_t i = 0; i < temps.size(); ++i) {
                if (temps[i] < stepMin) {
                    stepMin = temps[i];
                    stepMinNode = i + 1;
                }
                if (temps[i] > stepMax) {
                    stepMax = temps[i];
                    stepMaxNode = i + 1;
                }
            }
            
            std::cout << std::setw(7) << std::fixed << std::setprecision(1) << time << " | "
                     << std::setw(8) << std::setprecision(4) << stepMin << " (Node " 
                     << std::setw(2) << stepMinNode << ") | "
                     << std::setw(8) << std::setprecision(4) << stepMax << " (Node " 
                     << std::setw(2) << stepMaxNode << ") |\n";
        }
        std::cout << "\n" << std::string(60, '=') << "\n\n";
    }
    
    if (config.shouldPrintSteadyState()) {
        std::cout << "STEADY-STATE ANALYSIS\n\n";
        std::cout << "Solving [H+Hbc]{T} = {P}\n\n";
        
        auto steadyTemperatures = eqSystem.solve();
        eqSystem.printSolution(steadyTemperatures);
        
        std::cout << "\n" << std::string(60, '=') << "\n\n";
    }
    
    std::cout << "Simulation complete.\n";
    
    return 0;
}
