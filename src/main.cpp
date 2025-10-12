#include "Grid.h"
#include "GridSelector.h"
#include <iostream>
#include <string>

int main() {
    std::cout << "MES Grid Reader" << std::endl;
    std::cout << "===============" << std::endl;
    
    // Create grid selector
    GridSelector gridSelector("grids/");
    
    if (!gridSelector.hasGridFiles()) {
        std::cout << "No grid files found. Please add grid files to the grids/ directory." << std::endl;
        return 1;
    }
    
    while (true) {
        // Use GridSelector to get filename
        std::string filename = gridSelector.selectGrid();
        
        if (filename.empty()) {
            std::cout << "Goodbye!" << std::endl;
            break;
        }
        
        // Create and load grid
        Grid grid;
        if (grid.loadFromFile(filename)) {
            // Display basic information
            grid.printSummary();
            grid.printGlobalData();
            
            // Show node ids and coordinates
            grid.printNodesIdAndCoordinates();
            
            // Additional information menu
            char showMore;
            std::cout << "\nShow more details? (y/n): ";
            std::cin >> showMore;
            
            if (showMore == 'y' || showMore == 'Y') {
                grid.printElements();
                grid.printBoundaryConditions();
            }
        } else {
            std::cout << "Failed to load grid file: " << filename << std::endl;
        }
        
        char continueChoice;
        std::cout << "\nLoad another grid? (y/n): ";
        std::cin >> continueChoice;
        
        if (continueChoice != 'y' && continueChoice != 'Y') {
            break;
        }
        
        // Refresh file list in case new files were added
        gridSelector.refreshFileList();
    }
    
    return 0;
}
