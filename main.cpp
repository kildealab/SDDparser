#include <iostream>

#include "SDDparser.h"


int main()
{
    SDDparser parser;


    if(!parser.load("SDDOutput_1.txt"))		//INPUT SDD FILE PATH HERE
    {
        std::cout << "Failed to load SDD file\n";
        return 1;
    }


    const Header& header = parser.getHeader();


    std::cout << "====================\n";
    std::cout << "Header\n";
    std::cout << "====================\n";

// Instead of printing out redundant info, print out meaning of different
// numerical fields and summarize the SDD header. Later repeat for the 
// SDD exposure data entries.

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


// Print out summary of exposure entry data, figure out what is important to
// A potential user of the SDDparser
    for(const auto& exposure : parser.getExposures())
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


    return 0;
}
