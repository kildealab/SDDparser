#include <iostream>
#include <fstream>
#include <filesystem>

#include "SDDparser.h"
#include "Karyogram.h"

int main(int argc, char* argv[]) 					// Variables in main() brackets allow for the second command line argument to be the SDD file path
{

    if(argc != 2)							// This if statement allows command line arguments for the SDD file path 
    {
	std::cerr << "Usage: " << argv[0] << " <SDD file>\n";
	return 1;
    }

    std::string argument = argv[1];    					// check if -h or --help argument was passed instead of SDDfileName
    if(argument == "-h" || argument == "--help")			// Provide user with help on how to use the SDDparser in the command line.
    {
        std::cout
            << "SDDparser\n\n"
            << "Usage:\n"
            << "  " << argv[0] << " <SDD file>\n\n"
            << "Example:\n"
            << "  " << argv[0] << " SDDOutput_1.txt\n";
        return 0;
    }

    std::string filename = argv[1]; 					// SDDfile is the second argument of the command

    SDDparser parser; 							// Instantiate SDDparser object called parser

    if(!parser.load(filename))						// Check if SDD input file loads correctly
    {
        std::cout << "Failed to load SDD file " << filename << "\n";
        return 1;
    }

    std::filesystem::path inputPath(argv[1]);				// Create Summary Output file to the same path as the input SDD file

    std::filesystem::path summaryPath = inputPath.parent_path() / 	// Create the path to the summary file as the same path as the input SDD file
	(inputPath.stem().string() + "_summary.txt");

    std::ofstream summaryFile(summaryPath);				// Output to desired summary file path

    if(!summaryFile)							// Make sure summary file created successfully.
    {
        std::cerr
            << "Failed to create summary file\n";

        return 1;
    }

    parser.printSummary(summaryFile); 					//Returns Summary of SDD header fields in a summary file at the desired path

    summaryFile.close();

    std::cout << "Summary written to: " << summaryPath << "\n";		// Let user know summary file was created successfully

//// Karyogram drawing
    const Header& header = parser.getHeader();

    Karyogram karyogram;

    if (!karyogram.generateKaryogram(
            header.chromosome_sizes,
            header.cell_cycle_phase,
	    "karyogram.png"))
    {
        std::cerr << "Failed to generate karyogram.\n";
        return 1;
    }

    std::cout << "Karyogram generated successfully.\n";

////

    return 0;
}

