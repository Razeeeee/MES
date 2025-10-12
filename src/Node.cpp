#include "Node.h"
#include <iomanip>

Node::Node() : id(0), x(0.0), y(0.0) {}

Node::Node(int id, double x, double y) : id(id), x(x), y(y) {}

void Node::print() const {
    std::cout << "Node " << std::setw(4) << id 
              << ": (" << std::setw(12) << std::fixed << std::setprecision(6) << x 
              << ", " << std::setw(12) << std::fixed << std::setprecision(6) << y << ")" << std::endl;
}

std::ostream& operator<<(std::ostream& os, const Node& node) {
    os << "Node " << std::setw(4) << node.id 
       << ": (" << std::setw(12) << std::fixed << std::setprecision(6) << node.x 
       << ", " << std::setw(12) << std::fixed << std::setprecision(6) << node.y << ")";
    return os;
}
