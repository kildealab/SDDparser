#ifndef SDD_PARSER_H
#define SDD_PARSER_H

#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "SDDtypes.h"

class SDDparser
{

public:

    bool load(const std::string& filename);
    
    const Header& getHeader() const;

    const std::vector<Exposure>& getExposures() const;

private:

    Header header;

    bool parseHeader(std::ifstream& file);

    void parseDamage(std::ifstream& file);

    std::vector<Exposure> exposures;


};

#endif
