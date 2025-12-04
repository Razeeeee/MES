#include "GridSelector.h"
#include <iostream>
#include <algorithm>

GridSelector::GridSelector(const std::string& gridDir) 
    : gridDirectory(gridDir) {
    refreshFileList();
}

void GridSelector::scanGridDirectory() {
    availableGridFiles.clear();
    
    try {
        if (!std::filesystem::exists(gridDirectory)) {
            std::cerr << "Warning: Directory '" << gridDirectory << "' not found\n";
            return;
        }
        
        if (!std::filesystem::is_directory(gridDirectory)) {
            std::cerr << "Warning: '" << gridDirectory << "' is not a directory\n";
            return;
        }
        
        for (const auto& entry : std::filesystem::directory_iterator(gridDirectory)) {
            if (isValidGridFile(entry.path())) {
                availableGridFiles.push_back(entry.path().filename().string());
            }
        }
        
        // Sort files alphabetically
        std::sort(availableGridFiles.begin(), availableGridFiles.end());
        
    } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << "Error: Cannot access directory: " << ex.what() << "\n";
    }
}

bool GridSelector::isValidGridFile(const std::filesystem::path& filePath) const {
    return filePath.extension() == ".txt" && 
           std::filesystem::is_regular_file(filePath);
}

void GridSelector::displayMenu() const {
    std::cout << "\nSelect grid file:\n";
    
    for (size_t i = 0; i < availableGridFiles.size(); ++i) {
        std::cout << (i + 1) << ". " << availableGridFiles[i] << "\n";
    }
    
    std::cout << "0. Exit\n";
}

std::string GridSelector::selectGrid() {
    while (true) {
        displayMenu();
        
        int choice;
        std::cout << "\nChoice: ";
        
        if (!(std::cin >> choice)) {
            std::cout << "Invalid input. Enter a number.\n";
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            continue;
        }
        
        if (choice == 0) {
            return ""; // Empty string indicates exit
        }
        
        std::string filename = getGridPath(choice);
        if (filename.empty()) {
            std::cout << "Invalid choice. Try again.\n";
            continue;
        }
        
        return filename;
    }
}

std::string GridSelector::getGridPath(int choice) const {
    if (choice < 1 || choice > static_cast<int>(availableGridFiles.size())) {
        return "";
    }
    
    return gridDirectory + "/" + availableGridFiles[choice - 1];
}

bool GridSelector::hasGridFiles() const {
    return !availableGridFiles.empty();
}

void GridSelector::refreshFileList() {
    scanGridDirectory();
}
