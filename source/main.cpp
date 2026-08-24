#include <iostream>
#include <fstream>
#include <filesystem>

#include "SDDparser.h"
#include "Karyogram.h"

int main(int argc, char* argv[]) 					// Variables in main() brackets allow for the second command line argument to be the SDD file path or -h|--help for usage help. Third argument is an option for plotting the corresponding Karyogram if desired using "--karyogram". The fourth argument specifies either 'human' genome or 'other' genome.
{

    // ---------------------------------------
    // Check Command-line arguments
    // ---------------------------------------

    if(argc < 2)							// Must have at least two command-line arguments, the program and the SDD file path.
    {
	std::cerr << "Usage: ./SDDparser <SDD file> [--karyogram human|other]\n";
	return 1;
    }

    std::string argument = argv[1];					// If second command-line argument is -h or --help, return usage examples.
    if (argument == "-h" || argument == "--help")
    {
        std::cout
            << "SDDparser\n\n"
            << "Usage:\n"
            << "  ./SDDparser <SDD file>\n"
            << "  ./SDDparser <SDD file> --karyogram human\n"
            << "  ./SDDparser <SDD file> --karyogram other\n\n"
            << "Examples:\n"
            << "  ./SDDparser SDDOutput_1.txt\n"
            << "  ./SDDparser SDDOutput_1.txt --karyogram human\n"
            << "  ./SDDparser SDDOutput_1.txt --karyogram other\n";

        return 0;
    }

    std::string filename = argv[1]; 					// SDD file path is the second argument of the command-line

    // --------------------------------------
    // Karyogram options
    // --------------------------------------

    bool drawKaryogram = false;
    bool humanGenome = false;


    std::string genomeType = argv[3];					// The fourth argument is the genomeType = 'human'|'other'

    if (argc >= 3)							// If more at least 3 arguments are passed, must be '--karyogram'
    {
        if (std::string(argv[2]) == "--karyogram")
        {
            drawKaryogram = true;

            if (argc < 4)						// If '--karyogram' specified, user must specify either 'human' genome or 'other' genome for centromere plotting.
            {
                std::cerr << "ERROR: --karyogram requires "
                          << "either 'human' or 'other'.\n";
                return 1;
            }

            if (genomeType == "human")
            {
                humanGenome = true;
            }
            else if (genomeType == "other")
            {
                humanGenome = false;
            }
            else							// If fourth argument not 'human'|'other', return error.
            {
                std::cerr << "ERROR: Genome type must be "
                          << "'human' or 'other'.\n";
                return 1;
            }
	    if (argc > 4)						// Cannot accept more than 4 command-line arguments.
    	    {
		std::cerr << "Error: Too many command-line arguments.\n";
		return 1;
	    }
        }
        else								// Return error if unknown third command-line argument encountered.
        {
            std::cerr << "ERROR: Unknown command-line option: "
                      << argv[2] << '\n';
            return 1;
        }
    }

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

