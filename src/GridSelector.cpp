#include "GridSelector.h"
#include <iostream>
#include <algorithm>

GridSelector::GridSelector(const std::string& gridDir) 
    : gridDirectory(gridDir) {
    initializeDefaultFiles();
    refreshFileList();
}

void GridSelector::initializeDefaultFiles() {
    defaultGridFiles = {
        "Test1_4_4.txt",
        "Test2_4_4_MixGrid.txt", 
        "Test3_31_31_kwadrat.txt"
    };
}

void GridSelector::scanGridDirectory() {
    availableGridFiles.clear();
    
    try {
        if (std::filesystem::exists(gridDirectory) && 
            std::filesystem::is_directory(gridDirectory)) {
            
            for (const auto& entry : std::filesystem::directory_iterator(gridDirectory)) {
                if (isValidGridFile(entry.path())) {
                    availableGridFiles.push_back(entry.path().filename().string());
                }
            }
            
            // Sort files alphabetically
            std::sort(availableGridFiles.begin(), availableGridFiles.end());
        }
    } catch (const std::filesystem::filesystem_error& ex) {
        std::cout << "Warning: Could not access grid directory: " << ex.what() << std::endl;
        // Fall back to default files if directory is not accessible
        availableGridFiles = defaultGridFiles;
    }
    
    // If no files found in directory, use defaults
    if (availableGridFiles.empty()) {
        availableGridFiles = defaultGridFiles;
    }
}

bool GridSelector::isValidGridFile(const std::filesystem::path& filePath) const {
    return filePath.extension() == ".txt" && 
           std::filesystem::is_regular_file(filePath);
}

void GridSelector::displayMenu() const {
    std::cout << "\nSelect grid file:" << std::endl;
    
    for (size_t i = 0; i < availableGridFiles.size(); ++i) {
        std::cout << (i + 1) << ". " << availableGridFiles[i] << std::endl;
    }
    
    std::cout << "0. Exit" << std::endl;
}

std::string GridSelector::selectGrid() {
    while (true) {
        displayMenu();
        
        int choice;
        std::cout << "\nEnter your choice: ";
        
        if (!(std::cin >> choice)) {
            std::cout << "Invalid input. Please enter a number." << std::endl;
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            continue;
        }
        
        if (choice == 0) {
            return ""; // Empty string indicates exit
        }
        
        std::string filename = getGridPath(choice);
        if (filename.empty()) {
            std::cout << "Invalid choice. Please try again." << std::endl;
            continue;
        }
        
        return filename;
    }
}

std::string GridSelector::getGridPath(int choice) const {
    if (choice < 1 || choice > static_cast<int>(availableGridFiles.size())) {
        return "";
    }
    
    return gridDirectory + availableGridFiles[choice - 1];
}

bool GridSelector::hasGridFiles() const {
    return !availableGridFiles.empty();
}

void GridSelector::refreshFileList() {
    scanGridDirectory();
}
