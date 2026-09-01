#include <iostream>
#include <fstream>
#include <filesystem>

#include "SDDparser.h"
#include "SDRparser.h"
#include "Karyogram.h"


namespace
{
    void printUsage()
    {
        std::cout
            << "SDDparser\n\n"
            << "Usage:\n"
            << "  ./SDDparser <SDD file>\n"
            << "  ./SDDparser <SDD file> --karyogram human|other\n"
            << "  ./SDDparser -sdr <SDR file>\n"
            << "  ./SDDparser -sdr <SDR file> --karyogram human|other\n\n"
            << "Examples:\n"
            << "  ./SDDparser SDDOutput_1.txt\n"
            << "  ./SDDparser SDDOutput_1.txt --karyogram human\n"
            << "  ./SDDparser -sdr SDROutput_1.txt\n"
            << "  ./SDDparser -sdr SDROutput_1.txt --karyogram human\n";
    }
    
    // Parses an optional "--karyogram human|other" pair starting at
    // argv[startIndex]. If nothing follows the preceding arguments,
    // drawKaryogram is left false and this returns true (nothing to do).
    // Prints an error and returns false on any malformed input.
    bool parseKaryogramOption(
	int argc, 						// Number of command line arguments
	char* argv[], 						// The command line argument at a particular index
	int startIndex, 					// Variable to track the different possible commands
	bool& drawKaryogram, 					// Check if draw karyogram was desired by the user
	std::string& genomeType)				// Check if genome specified was human or other
    {
	drawKaryogram = false;
        genomeType.clear();

        if (argc <= startIndex)
        {
            return true;
        }

        if (std::string(argv[startIndex]) != "--karyogram")
        {
            std::cerr << "ERROR: Unknown command-line option: "
                      << argv[startIndex] << '\n';
            return false;
        }

        drawKaryogram = true;

	if (argc <= startIndex + 1)
        {
            std::cerr << "ERROR: --karyogram requires "
                      << "either 'human' or 'other'.\n";
            return false;
        }

        genomeType = argv[startIndex + 1];

        if (genomeType != "human" && genomeType != "other")
        {
            std::cerr << "ERROR: Genome type must be "
                      << "'human' or 'other'.\n";
            return false;
        }

        if (argc > startIndex + 2)
        {
            std::cerr << "ERROR: Too many command-line arguments.\n";
            return false;
        }

	return true;
    }
}










int main(int argc, char* argv[]) 					// Variables in main() brackets allow for the second command line argument to be the SDD file path or -h|--help for usage help. Third argument is an option for plotting the corresponding Karyogram if desired using "--karyogram". The fourth argument specifies either 'human' genome or 'other' genome.
{

    // ---------------------------------------
    // Check Command-line arguments
    // ---------------------------------------

    if(argc < 2)							// Must have at least two command-line arguments, the program and the SDD file path.
    {
	std::cerr << "Usage: ./SDDparser <SDD file> [--karyogram human|other]\n"
		   << "       ./SDDparser -sdr <SDR file> [--karyogram human|other]\n";
	return 1;
    }

    std::string firstArg = argv[1];					// If second command-line argument is -h or --help, return usage examples.
    if (firstArg == "-h" || firstArg == "--help")
    {
        printUsage();
        return 0;
    }

    // --------------------------------------
    // Check for -sdr file option
    // --------------------------------------

    if (firstArg == "-sdr")
    {
        if (argc < 3)
        {
            std::cerr << "ERROR: -sdr requires an SDR file path.\n";
            return 1;
        }

        std::string filename = argv[2];

        bool drawKaryogram = false;
        std::string genomeType;

	if (!parseKaryogramOption(argc, argv, 3, drawKaryogram, genomeType))
        {
            return 1;
        }

	bool humanGenome = (genomeType == "human");

        SDRparser parser;

        if (!parser.parseFile(filename))
        {
            std::cout << "Failed to load SDR file " << filename << "\n";
            return 1;
        }

	std::filesystem::path inputPath(filename);
        std::filesystem::path summaryPath = inputPath.parent_path() /
            (inputPath.stem().string() + "_summary.txt");

        if (!parser.writeSummary(summaryPath.string()))
        {
            std::cerr << "Failed to create SDR summary file\n";
            return 1;
        }

        std::cout << "Summary written to: " << summaryPath << "\n";

        if (drawKaryogram)
        {
	    Karyogram karyogram;
/*
	    std::filesystem::path karyogramPath =
            inputPath.parent_path() /
            (inputPath.stem().string() + "_karyogram_" + genomeType + ".png");
*/

	    bool hasFailed = false;

	    for (const SDRsubHeader& subHeader : parser.getSubHeaders())
	    {
		std::filesystem::path cellKaryogramPath =
            		inputPath.parent_path() /
            		(inputPath.stem().string() + "_cell" + std::to_string(subHeader.cellID) +
             		"_karyogram_" + genomeType + ".png");

		if (!karyogram.generateSDRkaryogram(parser.getMasterHeader(), subHeader, humanGenome, cellKaryogramPath.string()))
		{
		    std::cerr << "Failed to generate karyogram for cell " << subHeader.cellID << ".\n";
            	    hasFailed = true;
            	    continue;
		}

	        std::cout << "Karyogram generated successfully: " << cellKaryogramPath << "\n";

	    }

	    if (hasFailed)
	    {
		return 1;
	    }

        }

        return 0;
    }


    // --------------------------------------
    // SDD mode (existing behavior)
    // --------------------------------------

    std::string filename = firstArg;

    bool drawKaryogram = false;
    std::string genomeType;

    if (!parseKaryogramOption(argc, argv, 2, drawKaryogram, genomeType))
    {
	return 1;
    }

    bool humanGenome = (genomeType == "human");

    // -------------------------------------
    // Load SDD file
    // -------------------------------------

    SDDparser parser; 							// Instantiate SDDparser object called parser

    if(!parser.load(filename))						// Check if SDD input file loads correctly
    {
        std::cout << "Failed to load SDD file " << filename << "\n";
        return 1;
    }

    // -------------------------------------
    // Create summary file
    // -------------------------------------

    std::filesystem::path inputPath(filename);				// Create Summary Output file to the same path as the input SDD file

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

    // ---------------------------------------------
    // Karyogram drawing
    // ---------------------------------------------

    if (drawKaryogram)							// Check if user passed '--karyogram', meaning drawKaryogram == true. 
    {
        const Header& header = parser.getHeader();			// Get header to pass the chromosome sizes and cell cycle phase to the generate Karyogram function.

        Karyogram karyogram;						// Instantiate the Karyogram object called karyogram.

	std::filesystem::path karyogramPath =				// Create karyogram output filename using the input filename and generic suffix "_karyogram.png"
        inputPath.parent_path() /
        (inputPath.stem().string() + "_karyogram_" +  genomeType + ".png");


        if (!karyogram.generateKaryogram(				// Ensure Karyogram could be generated, if there are insufficient data in the SDD file.
                header.chromosome_sizes,
                header.cell_cycle_phase,
		header.dose_or_fluence,
		header.incident_particles,
	        parser.getExposures(),
		humanGenome,
	        karyogramPath.string()))
         {
            std::cerr << "Failed to generate karyogram.\n";
            return 1;
        }

        std::cout << "Karyogram generated successfully: " << karyogramPath << "\n"; // If SDD has sufficient data, exit with success message and generate the karyogram image.
    }


    return 0;
}

