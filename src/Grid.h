#pragma once
#include "Node.h"
#include "Element.h"
#include "GlobalData.h"
#include "EquationSystem.h"
#include <vector>
#include <string>
#include <set>

class Grid {
private:
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::set<int> boundaryConditions;
    GlobalData globalData;
    std::string filename;

    void parseGlobalData(const std::string& line);
    void parseNode(const std::string& line);
    void parseElement(const std::string& line);
    void parseBoundaryConditions(const std::string& line);
    
    std::string trim(const std::string& str);
    std::vector<std::string> split(const std::string& str, char delimiter);

public:
    Grid();
    Grid(const std::string& filename);
    
    bool loadFromFile(const std::string& filename);
    
    const std::vector<Node>& getNodes() const { return nodes; }
    const std::vector<Element>& getElements() const { return elements; }
    const std::set<int>& getBoundaryConditions() const { return boundaryConditions; }
    const GlobalData& getGlobalData() const { return globalData; }
    const std::string& getFilename() const { return filename; }
    
    Node* findNodeById(int id);
    Element* findElementById(int id);
    
    EquationSystem assembleGlobalEquationSystem(
        int numGaussPointsH = 2,
        int numGaussPointsC = 2, 
        int numGaussPointsBoundary = 2) const;
    
    std::vector<bool> getElementBoundaryEdges(const Element& element) const;
    
    void printGlobalData() const;
    void printNodes() const;
    void printElements() const;
    void printBoundaryConditions() const;
    void printSummary() const;
    void printNodesIdAndCoordinates() const;
};
