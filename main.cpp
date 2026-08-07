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


    std::cout << "Version: "
              << header.version
              << "\n";


    std::cout << "Software: "
              << header.software
              << "\n";


    std::cout << "Dose: "
              << header.dose
              << "\n";


    std::cout << "\nNumber of exposures: "
              << parser.getExposures().size()
              << "\n";


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
