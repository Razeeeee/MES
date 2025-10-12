#pragma once
#include "Node.h"
#include "Element.h"
#include "GlobalData.h"
#include <vector>
#include <string>
#include <set>

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
    std::set<int> boundaryConditions;
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
    const GlobalData& getGlobalData() const { return globalData; }
    const std::string& getFilename() const { return filename; }
    
    // Utility methods
    Node* findNodeById(int id);
    Element* findElementById(int id);
    
    // Display methods
    void printGlobalData() const;
    void printNodes() const;
    void printElements() const;
    void printBoundaryConditions() const;
    void printSummary() const;
    
    // Display node IDs and coordinates as requested
    void printNodesIdAndCoordinates() const;
};
