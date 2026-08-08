#include <iostream>
#include <fstream>

#include "SDDparser.h"


int main(int argc, char* argv[]) // Remove inside main if hardcoding filepath
{
//This if statement allows command line arguments for the SDD file path 
    if(argc != 2)//
    {//
	std::cerr << "Usage: " << argv[0] << " <SDD file>\n";//
	return 1;//
    }//

    std::string argument = argv[1];    
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

    SDDparser parser;

    if(!parser.load(filename))	//If command line arguments unwanted, remove above if statement and hardcode INPUT SDD FILE PATH HERE
    {
        std::cout << "Failed to load SDD file " << filename << "\n";
        return 1;
    }

    std::ofstream summaryFile(filename + "_header_summary.txt");

    if(!summaryFile)
    {
        std::cerr
            << "Failed to create summary file\n";

        return 1;
    }

    parser.printHeaderSummary(summaryFile); //Returns Summary of SDD header fields in a summary file

    summaryFile.close();
//    const Header& header = parser.getHeader();


//    std::cout << "====================\n";
//    std::cout << "Header\n";
//    std::cout << "====================\n";

// Instead of printing out redundant info, print out meaning of different
// numerical fields and summarize the SDD header. Later repeat for the 
// SDD exposure data entries.
/*
    std::cout << "SDD Version: "
              << header.sdd_version
              << "\n";


    std::cout << "Software: "
              << header.software
              << "\n";


    std::cout << "Dose or fluence: ";
    for (double val : header.dose_or_fluence) {
        std::cout << val << " ";
    }
    std::cout << "\n";


    std::cout << "\nNumber of exposures: "
              << parser.getExposures().size()
              << "\n";
*/

// Print out summary of exposure entry data, figure out what is important to
// A potential user of the SDDparser
/*    for(const auto& exposure : parser.getExposures())
    {
        std::cout
            << "\nExposure "
            << exposure.exposureID
            << " contains "
            << exposure.damages.size()
            << " entries\n";


        for(const auto& damage : exposure.damages)
        {
            std::cout
                << "  Event "
                << damage.classification.eventID
                << "\n";

            std::cout
                << "  Position: "
                << damage.position.x << ", "
                << damage.position.y << ", "
                << damage.position.z
                << "\n";
        }
    }

*/
    return 0;
}

