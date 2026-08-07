#ifndef SDD_PARSER_H
#define SDD_PARSER_H

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <ostream>
#include <iostream>

#include "SDDtypes.h"

class SDDparser
{

public:

    bool load(const std::string& filename);
    
    const Header& getHeader() const;

    const std::vector<Exposure>& getExposures() const;

    void printHeaderSummary(std::ostream& out = std::cout) const;

    void printDataEntrySummary(std::ostream& out = std::cout) const;

private:

    Header header;

    bool parseHeader(std::ifstream& file);

    void parseDamage(std::ifstream& file);

    std::vector<Exposure> exposures;


};

#endif
