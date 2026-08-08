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
	 std::cerr << "ERROR: SDD header parsing failed: "
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

// List all valid SDD header field names.
const std::unordered_set<std::string> validSDDHeaderFields =
{
    "sdd version", 			// Required
    "software",  			// Not required
    "author",				// Required
    "simulation details",   		// Not Required
    "source",  				// Required
    "source type", 			// Required
    "incident particles",  		// Not Required
    "mean particle energy",		// Not required
    "energy distribution",		// Not required
    "particle fraction",		// Not required
    "dose or fluence",			// Not required
    "dose rate",			// Not required
    "irradiation target",		// Not Required
    "volumes",				// Required
    "chromosome sizes",			// Required
    "dna density",			// Required
    "cell cycle phase",			// Required
    "dna structure",			// Required
    "in vitro / in vivo",		// Not required
    "proliferation status",		// Not required
    "microenvironment",			// Not required
    "damage definition",		// Required
    "time",				// Not required
    "damage and primary count",		// Required
    "data entries",			// Required
    "additional information"		// Not required
};

bool SDDparser::parseHeader(std::ifstream& file)
{
    std::string line;


    while(std::getline(file,line))
    {

        line = trim(line);


        if(line.empty())
        {
	    continue;
	}
        
	// Remove final semicolon
        if(line.back()==';')
            line.pop_back();

        if(normalizeHeaderKey(line) == "***endofheader***")
        {    std::cout << "SDD header parsed successfully.\n";
		
	     return true;
	}

        // Split header entry by commas
        std::vector<std::string> tokens = split(line,',');

        if(tokens.empty())
            continue;

        std::string original_key = tokens[0]; //Before removing whitespaces and capitalizations
	std::string key = normalizeHeaderKey(original_key); // Headers are now lower case and have no leading or trailing whitespaces

	// Printing out current SDD header field being read for DEBUGGING
	std::cout << "Reading header field: [" << original_key << "]\n";

	// Check if SDD header field is supported
	if(validSDDHeaderFields.find(key) == validSDDHeaderFields.end())
	{
    	    std::cerr
            << "ERROR: Unknown SDD header field: " << original_key << std::endl;
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

	// Fields with different numbering data types are converted to floats.
	// Fields with strings and ints/floats are kept as strings.

        if(key=="sdd version")
        {
            header.sdd_version = values[0];
        }

        else if(key=="software")
        {
            header.software = values[0];
        }

        else if(key=="author")
        {
            header.author = values[0];
        }

	else if(key=="simulation details")
	{
	    header.simulation_details = values[0]; 
	}

	else if(key=="source")
	{
	    header.source = values[0];
	}

	else if(key=="source type")
	{
	    header.source_type = std::stoi(values[0]);
	}

	else if(key=="incident particles")
	{
 	    header.incident_particles = parseIntList(values); 
	}

	else if(key=="mean particle energy")
	{
	    header.mean_particle_energy = parseDoubleList(values);
	}

	else if(key=="energy distribution")
	{
	    header.energy_distribution = values[0];
	}

	else if(key=="particle fraction")
	{
	    header.particle_fraction = parseDoubleList(values);
	}

        else if(key=="dose or fluence")
        {
            header.dose_or_fluence = parseDoubleList(values);
	}

        else if(key=="dose rate")
        {
            header.dose_rate = std::stod(values[0]);
        }

	else if(key=="irradiation target")
	{
	    header.irradiation_target = values[0];
	}

	else if(key=="volumes")
	{
	    header.volumes = parseDoubleList(values);
	}

        else if(key=="chromosome sizes")
        {
            header.chromosome_sizes = parseDoubleList(values);
        }

	else if(key=="dna density")
	{
            header.DNA_density = std::stod(values[0]);
	}

	else if(key=="cell cycle phase")
	{
	    header.cell_cycle_phase = parseDoubleList(values);
	}	

	else if(key=="dna structure")
	{
	    header.DNA_structure = parseIntList(values);
	}

	else if(key=="in vitro / in vivo")
	{
	    header.in_vitro_or_in_vivo = std::stoi(values[0]);
	}

	else if(key=="proliferation status")
	{
	    header.proliferation_status = std::stoi(values[0]);
	}
	
	else if(key=="microenvironment")
	{
	    header.microenvironment = parseDoubleList(values);
	}
	
	else if(key=="damage definition")
	{
	    header.damage_definition = parseDoubleList(values);
	}

	else if(key=="time")
	{
	    header.time = std::stod(values[0]);
	}

	else if(key=="damage and primary count")
	{
	    header.damage_and_primary_count = parseIntList(values);
	}

	else if(key=="data entries")
	{
	    header.data_entries = parseIntList(values);
	}

	else if(key=="additional information")
	{
	    header.additional_information = values[0];
	}

        else
        {
            // Unknown SDD field
            std::cerr << "Unknown SDD header field entry = " << key << "\n";

        }

    }
	// If ***EndOfHeader*** not found, exit with error
	std::cerr << "ERROR: End of SDD header not found.\n";
	
	return false;
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


void SDDparser::printHeaderSummary(std::ostream& out) const {
// MODIFY TO PRINT HEADER FIELDS
	out << "This is a test that SDDparser reads the SDD header ";
}

