#include <iostream>
#include <fstream>
#include <filesystem>

#include "SDDparser.h"


int main(int argc, char* argv[]) // Remove inside main if hardcoding filepath
{
////This if statement allows command line arguments for the SDD file path 
    if(argc != 2)
    {
	std::cerr << "Usage: " << argv[0] << " <SDD file>\n";
	return 1;
    }
////
    std::string argument = argv[1];    // check if -h or --help argument was passed instead of SDDfileName
    if(argument == "-h" || argument == "--help")
    {
        std::cout
            << "SDDparser\n\n"
            << "Usage:\n"
            << "  " << argv[0] << " <SDD file>\n\n"
            << "Example:\n"
            << "  " << argv[0] << " SDDOutput_1.txt\n";

        return 0;
    }

    std::string filename = argv[1]; // SDDfile is the second argument of the command

    SDDparser parser; 		// Instantiate SDDparser object called parser

// Check if SDD input file loads correctly
    if(!parser.load(filename))	//If command line arguments unwanted, remove above if statement and hardcode INPUT SDD FILE PATH in the place of "filename"
    {
        std::cout << "Failed to load SDD file " << filename << "\n";
        return 1;
    }

// Create Summary Output file to the same path as the input SDD file
    std::filesystem::path inputPath(argv[1]);

    std::filesystem::path summaryPath =
    inputPath.parent_path() /
    (inputPath.stem().string() + "_summary.txt");

    std::ofstream summaryFile(summaryPath);

    if(!summaryFile)
    {
        std::cerr
            << "Failed to create summary file\n";

        return 1;
    }

    parser.printSummary(summaryFile); //Returns Summary of SDD header fields in a summary file

    summaryFile.close();

    std::cout << "Summary written to: " << summaryPath << "\n";

    return 0;
}

