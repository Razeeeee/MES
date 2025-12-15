#include "Node.h"
#include <iomanip>

Node::Node() : id(0), x(0.0), y(0.0), hasBC(false) {}

Node::Node(int id, double x, double y, bool hasBC) : id(id), x(x), y(y), hasBC(hasBC) {}

void Node::print() const {
    std::cout << *this << std::endl;
}

std::ostream& operator<<(std::ostream& os, const Node& node) {
    os << "Node " << std::setw(4) << node.id 
       << ": (" << std::setw(12) << std::fixed << std::setprecision(6) << node.x 
       << ", " << std::setw(12) << std::fixed << std::setprecision(6) << node.y << ")";
    if (node.hasBC) os << " [BC]";
    return os;
}
