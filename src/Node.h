#pragma once
#include <iostream>

/**
 * @brief Represents a node in the finite element mesh
 * 
 * Contains node ID and 2D coordinates (x, y)
 */
class Node {
private:
    int id;
    double x, y;

public:
    Node();
    Node(int id, double x, double y);
    
    // Getters
    int getId() const noexcept { return id; }
    double getX() const noexcept { return x; }
    double getY() const noexcept { return y; }
    
    // Setters
    void setId(int id) noexcept { this->id = id; }
    void setX(double x) noexcept { this->x = x; }
    void setY(double y) noexcept { this->y = y; }
    
    // Display node information
    void print() const;
    
    // Friend operator for output stream
    friend std::ostream& operator<<(std::ostream& os, const Node& node);
};
