#include "GlobalData.h"
#include <iostream>
#include <iomanip>

GlobalData::GlobalData() 
    : simulationTime(0.0), simulationStepTime(0.0), conductivity(0.0),
      alfa(0.0), tot(0.0), initialTemp(0.0), density(0.0), specificHeat(0.0),
      nodesNumber(0), elementsNumber(0) {}

void GlobalData::print() const {
    std::cout << "\n=== GLOBAL DATA ===" << std::endl;
    std::cout << std::left;
    std::cout << std::setw(20) << "Simulation Time:" << simulationTime << std::endl;
    std::cout << std::setw(20) << "Step Time:" << simulationStepTime << std::endl;
    std::cout << std::setw(20) << "Conductivity:" << conductivity << std::endl;
    std::cout << std::setw(20) << "Alfa:" << alfa << std::endl;
    std::cout << std::setw(20) << "Tot:" << tot << std::endl;
    std::cout << std::setw(20) << "Initial Temp:" << initialTemp << std::endl;
    std::cout << std::setw(20) << "Density:" << density << std::endl;
    std::cout << std::setw(20) << "Specific Heat:" << specificHeat << std::endl;
    std::cout << std::setw(20) << "Nodes Number:" << nodesNumber << std::endl;
    std::cout << std::setw(20) << "Elements Number:" << elementsNumber << std::endl;
    std::cout << "===================" << std::endl;
}
