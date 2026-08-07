#include "SDDparser.h"
#include "SDDutilities.h"

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

bool SDDparser::load(const std::string& filename)
{
    std::ifstream file(filename);

    if(!file)
        { 
	 std::cerr << "ERROR: Cannot open file: "
		   << filename << std::endl ;

	 return false;
	}

    header = Header();
    exposures.clear();

    if(!parseHeader(file))
    {
    	std::cerr << "ERROR: SDD header parsing failed." << std::endl;
        return false;
    }

    parseDamage(file);

    return true;
}

// Add more header fields
const std::unordered_set<std::string> validSDDHeaderFields =
{
    "SDD Version",
    "Software",
    "Author",
    "Simulation details",
    "Irradiation target",
    "Source type",
    "Dose",
    "Dose rate",
    "Chromosome sizes",
    "Data entries"
};

bool SDDparser::parseHeader(std::ifstream& file)
{
    std::string line;


    while(std::getline(file,line))
    {

        line = trim(line);


        if(line.empty())
            continue;


        if(line == "***EndOfHeader***")
            break;

        // Remove final semicolon
        if(line.back()==';')
            line.pop_back();

        // Split header entry by commas
        std::vector<std::string> tokens = split(line,',');

        if(tokens.empty())
            continue;

        std::string key = tokens[0];


	// Printing out current SDD header field being read
	std::cout << "Reading header field: [" << key << "]\n";

	// Check if SDD header field is supported
	if(validSDDHeaderFields.find(key) == validSDDHeaderFields.end())
	{
    	    std::cerr
            << "ERROR: Unknown SDD header field: " << key << std::endl;
    	    return false;
	}

        // Remaining values
        std::vector<std::string> values(
            tokens.begin()+1,
            tokens.end()
        );

        //----------------------------
        // Map SDD fields
        //----------------------------


        if(key=="SDD Version")
        {
            header.version = values[0];
        }


        else if(key=="Software")
        {
            header.software = values[0];
        }


        else if(key=="Author")
        {
            header.author = values[0];
        }


        else if(key=="Dose")
        {
            header.dose =
                std::stod(values[0]);
        }


        else if(key=="Dose rate")
        {
            header.doseRate =
                std::stod(values[0]);
        }


        else if(key=="Chromosome sizes")
        {
            header.chromosomeSizes =
                parseDoubleList(values);
        }


        else if(key=="Data entries")
        {
            header.dataEntries =
                parseIntList(values);
        }


        else if(key=="Source type")
        {
            header.sourceType =
                values[0];
        }


        else if(key=="Irradiation target")
        {
            header.irradiationTarget =
                values[0];
        }


        else
        {
            // Unknown SDD field
            std::cerr << "Unknown SDD header field entry = " << key << "\n";
		// Ignore for now
        }

    }

	return true;
}


void SDDparser::parseDamage(std::ifstream& file)
{
    std::string line;

    while(std::getline(file,line))
    {
        if(line.empty())
            continue;


        DamageEntry damage;

        damage.rawLine = line;


        // TODO:
        // Parse SDD Data entries here
        // according to header.dataEntries


        Exposure exposure;

        exposure.exposureID =
            damage.classification.exposureID;

        exposure.damages.push_back(damage);

        exposures.push_back(exposure);
    }
}


const Header& SDDparser::getHeader() const
{
    return header;
}

const std::vector<Exposure>& SDDparser::getExposures() const 
{
    return exposures;
}
