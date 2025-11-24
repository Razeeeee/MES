#pragma once
#include "Node.h"
#include "Element.h"
#include "GlobalData.h"
#include "EquationSystem.h"
#include <vector>
#include <string>
#include <set>
#include <map>

/**
 * @brief Structure to hold boundary condition parameters for a node
 */
struct BoundaryCondition {
    double alfa;
    double tot;
    
    BoundaryCondition() : alfa(0.0), tot(0.0) {}
    BoundaryCondition(double a, double t) : alfa(a), tot(t) {}
};

/**
 * @brief Main grid class that manages the finite element mesh
 * 
 * Handles reading grid files, storing nodes and elements,
 * and providing access to the mesh data
 */
class Grid {
private:
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::set<int> boundaryConditions; // For backward compatibility - simple BC format
    std::map<int, BoundaryCondition> boundaryConditionsMap; // New format with node-specific parameters
    GlobalData globalData;
    std::string filename;

    // Helper methods for parsing
    void parseGlobalData(const std::string& line);
    void parseNode(const std::string& line);
    void parseElement(const std::string& line);
    void parseBoundaryConditions(const std::string& line);
    
    // String utility methods
    std::string trim(const std::string& str);
    std::vector<std::string> split(const std::string& str, char delimiter);

public:
    Grid();
    Grid(const std::string& filename);
    
    // File operations
    bool loadFromFile(const std::string& filename);
    
    // Getters
    const std::vector<Node>& getNodes() const { return nodes; }
    const std::vector<Element>& getElements() const { return elements; }
    const std::set<int>& getBoundaryConditions() const { return boundaryConditions; }
    const std::map<int, BoundaryCondition>& getBoundaryConditionsMap() const { return boundaryConditionsMap; }
    const GlobalData& getGlobalData() const { return globalData; }
    const std::string& getFilename() const { return filename; }
    
    // Get boundary condition parameters for a specific node
    // Returns true if node has BC, fills alfa and tot parameters
    bool getNodeBoundaryCondition(int nodeId, double& alfa, double& tot) const;
    
    // Utility methods
    Node* findNodeById(int id);
    Element* findElementById(int id);
    
    // Global matrix assembly
    // Calculate global H matrix by aggregating local element H matrices
    // Returns EquationSystem object containing the global H matrix [N x N]
    EquationSystem assembleGlobalEquationSystem(int numGaussPoints = 2) const;
    
    // Helper method to determine which edges of an element have boundary conditions
    // Returns vector of 4 bools: [bottom, right, top, left]
    std::vector<bool> getElementBoundaryEdges(const Element& element) const;
    
    // Display methods
    void printGlobalData() const;
    void printNodes() const;
    void printElements() const;
    void printBoundaryConditions() const;
    void printSummary() const;
    
    // Display node IDs and coordinates as requested
    void printNodesIdAndCoordinates() const;
};
