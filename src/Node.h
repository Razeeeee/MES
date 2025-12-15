#pragma once
#include <iostream>

class Node {
private:
    int id;
    double x, y;
    bool hasBC;

public:
    Node();
    Node(int id, double x, double y, bool hasBC = false);
    
    int getId() const noexcept { return id; }
    double getX() const noexcept { return x; }
    double getY() const noexcept { return y; }
    bool getBC() const noexcept { return hasBC; }
    
    void setId(int id) noexcept { this->id = id; }
    void setX(double x) noexcept { this->x = x; }
    void setY(double y) noexcept { this->y = y; }
    void setBC(bool bc) noexcept { this->hasBC = bc; }
    
    void print() const;
    
    friend std::ostream& operator<<(std::ostream& os, const Node& node);
};
