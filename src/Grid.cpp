#include "Grid.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <iomanip>

Grid::Grid() {}

Grid::Grid(const std::string& filename) {
    loadFromFile(filename);
}

bool Grid::loadFromFile(const std::string& filename) {
    this->filename = filename;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open " << filename << "\n";
        return false;
    }
    
    std::string line;
    std::string currentSection = "";
    
    // Clear existing data
    nodes.clear();
    elements.clear();
    boundaryConditions.clear();
    
    while (std::getline(file, line)) {
        line = trim(line);
        
        if (line.empty()) continue;
        
        // Check for section headers
        if (line.find("*Node") == 0) {
            currentSection = "NODE";
            continue;
        } else if (line.find("*Element") == 0) {
            currentSection = "ELEMENT";
            continue;
        } else if (line.find("*BC") == 0) {
            currentSection = "BC";
            continue;
        }
        
        // Parse data based on current section
        if (currentSection == "NODE") {
            parseNode(line);
        } else if (currentSection == "ELEMENT") {
            parseElement(line);
        } else if (currentSection == "BC") {
            parseBoundaryConditions(line);
        } else {
            // Parse global data (before any section)
            parseGlobalData(line);
        }
    }
    
    file.close();
    std::cout << "Loaded: " << filename << "\n";
    return true;
}

void Grid::parseGlobalData(const std::string& line) {
    std::istringstream iss(line);
    std::string key;
    double value;
    
    if (iss >> key >> value) {
        if (key == "SimulationTime") {
            globalData.setSimulationTime(value);
        } else if (key == "SimulationStepTime") {
            globalData.setSimulationStepTime(value);
        } else if (key == "Conductivity") {
            globalData.setConductivity(value);
        } else if (key == "Alfa") {
            globalData.setAlfa(value);
        } else if (key == "Tot") {
            globalData.setTot(value);
        } else if (key == "InitialTemp") {
            globalData.setInitialTemp(value);
        } else if (key == "Density") {
            globalData.setDensity(value);
        } else if (key == "SpecificHeat") {
            globalData.setSpecificHeat(value);
        } else if (key == "Nodes") {
            std::string number;
            int nodesNum;
            if (iss >> number >> nodesNum && number == "number") {
                globalData.setNodesNumber(nodesNum);
            }
        } else if (key == "Elements") {
            std::string number;
            int elementsNum;
            if (iss >> number >> elementsNum && number == "number") {
                globalData.setElementsNumber(elementsNum);
            }
        }
    }
}

void Grid::parseNode(const std::string& line) {
    // Remove leading/trailing whitespace and commas
    std::string cleanLine = line;
    std::replace(cleanLine.begin(), cleanLine.end(), ',', ' ');
    
    std::istringstream iss(cleanLine);
    int id;
    double x, y;
    
    if (iss >> id >> x >> y) {
        nodes.emplace_back(id, x, y);
    }
}

void Grid::parseElement(const std::string& line) {
    // Remove leading/trailing whitespace and commas
    std::string cleanLine = line;
    std::replace(cleanLine.begin(), cleanLine.end(), ',', ' ');
    
    std::istringstream iss(cleanLine);
    int id, nodeId;
    std::vector<int> nodeIds;
    
    if (iss >> id) {
        while (iss >> nodeId) {
            nodeIds.push_back(nodeId);
        }
        if (!nodeIds.empty()) {
            elements.emplace_back(id, nodeIds);
        }
    }
}

void Grid::parseBoundaryConditions(const std::string& line) {
    // Remove leading/trailing whitespace and commas
    std::string cleanLine = line;
    std::replace(cleanLine.begin(), cleanLine.end(), ',', ' ');
    
    std::istringstream iss(cleanLine);
    int nodeId;
    
    while (iss >> nodeId) {
        boundaryConditions.insert(nodeId);
    }
}

std::string Grid::trim(const std::string& str) {
    size_t first = str.find_first_not_of(' ');
    if (std::string::npos == first) {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

std::vector<std::string> Grid::split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(trim(token));
    }
    
    return tokens;
}

Node* Grid::findNodeById(int id) {
    for (auto& node : nodes) {
        if (node.getId() == id) {
            return &node;
        }
    }
    return nullptr;
}

Element* Grid::findElementById(int id) {
    for (auto& element : elements) {
        if (element.getId() == id) {
            return &element;
        }
    }
    return nullptr;
}

void Grid::printGlobalData() const {
    globalData.print();
}

void Grid::printNodes() const {
    std::cout << "\nNODES\n";
    for (const auto& node : nodes) {
        node.print();
    }
    std::cout << "Total: " << nodes.size() << "\n";
}

void Grid::printElements() const {
    std::cout << "\nELEMENTS\n";
    for (const auto& element : elements) {
        element.print();
    }
    std::cout << "Total: " << elements.size() << "\n";
}

void Grid::printBoundaryConditions() const {
    std::cout << "\nBOUNDARY CONDITIONS\n";
    std::cout << "Nodes: ";
    for (auto it = boundaryConditions.begin(); it != boundaryConditions.end(); ++it) {
        std::cout << *it;
        if (std::next(it) != boundaryConditions.end()) {
            std::cout << ", ";
        }
    }
    std::cout << "\n";
    std::cout << "Total: " << boundaryConditions.size() << " nodes\n";
}

void Grid::printSummary() const {
    std::cout << "\nGRID SUMMARY\n";
    std::cout << "File: " << filename << "\n";
    std::cout << "Nodes: " << nodes.size() << "\n";
    std::cout << "Elements: " << elements.size() << "\n";
    std::cout << "Boundary Nodes: " << boundaryConditions.size() << "\n";
    std::cout << std::string(40, '=') << "\n";
}

void Grid::printNodesIdAndCoordinates() const {
    std::cout << "\n=== NODE IDs AND COORDINATES ===" << std::endl;
    std::cout << std::left << std::setw(8) << "Node ID" 
              << std::setw(15) << "X Coordinate" 
              << std::setw(15) << "Y Coordinate" << std::endl;
    std::cout << std::string(38, '-') << std::endl;
    
    for (const auto& node : nodes) {
        std::cout << std::left << std::setw(8) << node.getId()
                  << std::setw(15) << std::fixed << std::setprecision(6) << node.getX()
                  << std::setw(15) << std::fixed << std::setprecision(6) << node.getY()
                  << std::endl;
    }
    std::cout << "=================================" << std::endl;
    std::cout << "Total nodes: " << nodes.size() << std::endl;
}

EquationSystem Grid::assembleGlobalEquationSystem(int numGaussPointsH, int numGaussPointsC, int numGaussPointsBoundary) const {
    (void)numGaussPointsC; // Unused - C matrix calculated per element
    // Get the number of nodes (size of global system matrices)
    int N = nodes.size();
    
    // Initialize equation system with size N
    EquationSystem eqSystem(N);
    
    // Get material and boundary properties from global data
    double conductivity = globalData.getConductivity();
    double alfa = globalData.getAlfa();
    double tot = globalData.getTot();
    double density = globalData.getDensity();
    double specificHeat = globalData.getSpecificHeat();
    
    // Assemble global matrices by iterating through all elements
    // FEM Assembly: Σ (local contributions) → global matrices
    for (const auto& element : elements) {
        // Get node IDs for this element (1-indexed in grid file)
        const auto& nodeIds = element.getNodeIds();
        
        // Extract node coordinates for this element [4 nodes]
        std::vector<double> nodeX(4), nodeY(4);
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            const Node& node = nodes[nodeIds[i] - 1]; // Convert to 0-indexed
            nodeX[i] = node.getX();
            nodeY[i] = node.getY();
        }
        
        // OPTIMIZED: Calculate H and C matrices together (single Jacobian calculation per integration point)
        auto [localH, localC] = element.calculateHAndCMatrices(
            nodeX, nodeY, conductivity, density, specificHeat, numGaussPointsH);
        
        // Determine which edges have convection boundary conditions
        auto boundaryEdges = getElementBoundaryEdges(element);
        
        // OPTIMIZED: Calculate Hbc and P together for boundary conditions
        auto [localHbc, localP] = element.calculateHbcAndPVector(
            nodeX, nodeY, alfa, tot, boundaryEdges, numGaussPointsBoundary);
        
        // Assemble local matrices into global system
        // Map from local element indices (0-3) to global node indices
        for (size_t i = 0; i < nodeIds.size(); ++i) {
            int globalI = nodeIds[i] - 1; // Convert to 0-indexed
            
            // Add local P contribution to global P vector
            eqSystem.addToPVector(globalI, localP[i]);
            
            for (size_t j = 0; j < nodeIds.size(); ++j) {
                int globalJ = nodeIds[j] - 1; // Convert to 0-indexed
                
                // Add local contributions to global matrices
                eqSystem.addToHMatrix(globalI, globalJ, localH[i][j]);
                eqSystem.addToHbcMatrix(globalI, globalJ, localHbc[i][j]);
                eqSystem.addToCMatrix(globalI, globalJ, localC[i][j]);
            }
        }
    }
    
    return eqSystem;
}

std::vector<bool> Grid::getElementBoundaryEdges(const Element& element) const {
    // Edge definitions for 4-node quadrilateral (DC2D4):
    // Edge 0 (bottom): nodes 1-2 (local 0-1)
    // Edge 1 (right):  nodes 2-3 (local 1-2)
    // Edge 2 (top):    nodes 3-4 (local 2-3)
    // Edge 3 (left):   nodes 4-1 (local 3-0)
    
    const auto& nodeIds = element.getNodeIds();
    std::vector<bool> boundaryEdges(4, false);
    
    // Check each edge
    // Edge 0 (bottom): nodes 0-1
    if (boundaryConditions.count(nodeIds[0]) && boundaryConditions.count(nodeIds[1])) {
        boundaryEdges[0] = true;
    }
    
    // Edge 1 (right): nodes 1-2
    if (boundaryConditions.count(nodeIds[1]) && boundaryConditions.count(nodeIds[2])) {
        boundaryEdges[1] = true;
    }
    
    // Edge 2 (top): nodes 2-3
    if (boundaryConditions.count(nodeIds[2]) && boundaryConditions.count(nodeIds[3])) {
        boundaryEdges[2] = true;
    }
    
    // Edge 3 (left): nodes 3-0
    if (boundaryConditions.count(nodeIds[3]) && boundaryConditions.count(nodeIds[0])) {
        boundaryEdges[3] = true;
    }
    
    return boundaryEdges;
}
