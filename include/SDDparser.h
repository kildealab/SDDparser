#ifndef SDD_PARSER_H
#define SDD_PARSER_H

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <ostream>
#include <iostream>
#include <fstream>

#include "SDDtypes.h"

class SDDparser
{

public:

    bool load(const std::string& filename);
    
    const Header& getHeader() const;

    const std::vector<Exposure>& getExposures() const;

// Function that will run printHeaderSummary and printDataEntrySummary
    void printSummary(std::ostream& out) const;

private:

    Header header;

    bool parseHeader(std::ifstream& file);

    void parseDamage(std::ifstream& file);

    std::vector<Exposure> exposures;

// Function to print the SDD header summary to an output file 
    void printHeaderSummary(std::ostream& out) const;

// std::cout to print summary to screen, remove to output to a file.
    void printExposureSummary(std::ostream& out = std::cout) const;

    std::map<int, std::size_t> exposureMap;

};

#endif
