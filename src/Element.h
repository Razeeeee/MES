#pragma once
#include <vector>
#include <iostream>

/**
 * @brief Represents a finite element in the mesh
 * 
 * Contains element ID and node IDs that form the element
 * Supports DC2D4 (4-node quadrilateral) elements
 */
class Element {
private:
    int id;
    std::vector<int> nodeIds;
    std::string type;

public:
    Element();
    Element(int id, const std::vector<int>& nodeIds, const std::string& type = "DC2D4");
    
    // Getters
    int getId() const { return id; }
    const std::vector<int>& getNodeIds() const { return nodeIds; }
    const std::string& getType() const { return type; }
    
    // Setters
    void setId(int id) { this->id = id; }
    void setNodeIds(const std::vector<int>& nodeIds) { this->nodeIds = nodeIds; }
    void setType(const std::string& type) { this->type = type; }
    
    // Display element information
    void print() const;
    
    // Friend operator for output stream
    friend std::ostream& operator<<(std::ostream& os, const Element& element);
};
