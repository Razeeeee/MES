#pragma once
#include <string>

/**
 * @brief Stores global simulation parameters
 * 
 * Contains all the global data parameters from the grid file
 * such as simulation time, material properties, etc.
 */
class GlobalData {
private:
    double simulationTime;
    double simulationStepTime;
    double conductivity;
    double alfa;
    double tot;
    double initialTemp;
    double density;
    double specificHeat;
    int nodesNumber;
    int elementsNumber;

public:
    GlobalData();
    
    // Getters
    double getSimulationTime() const { return simulationTime; }
    double getSimulationStepTime() const { return simulationStepTime; }
    double getConductivity() const { return conductivity; }
    double getAlfa() const { return alfa; }
    double getTot() const { return tot; }
    double getInitialTemp() const { return initialTemp; }
    double getDensity() const { return density; }
    double getSpecificHeat() const { return specificHeat; }
    int getNodesNumber() const { return nodesNumber; }
    int getElementsNumber() const { return elementsNumber; }
    
    // Setters
    void setSimulationTime(double value) { simulationTime = value; }
    void setSimulationStepTime(double value) { simulationStepTime = value; }
    void setConductivity(double value) { conductivity = value; }
    void setAlfa(double value) { alfa = value; }
    void setTot(double value) { tot = value; }
    void setInitialTemp(double value) { initialTemp = value; }
    void setDensity(double value) { density = value; }
    void setSpecificHeat(double value) { specificHeat = value; }
    void setNodesNumber(int value) { nodesNumber = value; }
    void setElementsNumber(int value) { elementsNumber = value; }
    
    // Display global data information
    void print() const;
};
