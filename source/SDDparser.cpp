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
    "Source",
    "Source type",
    "Incident particles",
    "Mean particle energy",
    "Energy distrbution",
    "Dose or fluence",
    "Irradiation target",
    "Volumes",
    "Dose rate",
    "Chromosome sizes",
    "DNA density",
    "Cell cycle phase",
    "DNA structure",
    "In vitro / in vivo",
    "Proliferation status",
    "Microenvironment",
    "Damage definition",
    "Time",
    "Damage and primary count",
    "Data entries",
    "Additional information"
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

	// Fields with different numbering data types are converted to floats.
	// Fields with strings and ints/floats are kept as strings.

        if(key=="SDD Version")
        {
            header.sdd_version = values[0];
        }

        else if(key=="Software")
        {
            header.software = values[0];
        }

        else if(key=="Author")
        {
            header.author = values[0];
        }

	else if(key=="Simulation details")
	{
	    header.simulation_details = values[0]; 
	}

	else if(key=="Source")
	{
	    header.source = values[0];
	}

	else if(key=="Source type")
	{
	    header.source_type = std::stoi(values[0]);
	}

	else if(key=="Incident particles")
	{
 	    header.incident_particles = parseIntList(values); 
	}

	else if(key=="Mean particle energy")
	{
	    header.mean_particle_energy = parseDoubleList(values);
	}

	else if(key=="Energy distribution")
	{
	    header.energy_distribution = values[0];
	}

	else if(key=="Particle fraction")
	{
	    header.particle_fraction = parseDoubleList(values);
	}

        else if(key=="Dose or fluence")
        {
            header.dose_or_fluence = parseDoubleList(values);
	}

        else if(key=="Dose rate")
        {
            header.dose_rate = std::stod(values[0]);
        }

	else if(key=="Irradiation target")
	{
	    header.irradiation_target = values[0];
	}

	else if(key=="Volumes")
	{
	    header.volumes = parseDoubleList(values);
	}

        else if(key=="Chromosome sizes")
        {
            header.chromosome_sizes = parseDoubleList(values);
        }

	else if(key=="DNA density")
	{
            header.DNA_density = std::stod(values[0]);
	}

	else if(key=="Cell cycle phase")
	{
	    header.cell_cycle_phase = parseDoubleList(values);
	}	

	else if(key=="DNA structure")
	{
	    header.DNA_structure = parseIntList(values);
	}

	else if(key=="In vitro / in vivo")
	{
	    header.in_vitro_or_in_vivo = std::stoi(values[0]);
	}

	else if(key=="Proliferation status")
	{
	    header.proliferation_status = std::stoi(values[0]);
	}
	
	else if(key=="Microenvironment")
	{
	    header.microenvironment = parseDoubleList(values);
	}
	
	else if(key=="Damage definition")
	{
	    header.damage_definition = parseDoubleList(values);
	}

	else if(key=="Time")
	{
	    header.time = std::stod(values[0]);
	}

	else if(key=="Damage and primary count")
	{
	    header.damage_and_primary_count = parseIntList(values);
	}

	else if(key=="Data entries")
	{
	    header.data_entries = parseIntList(values);
	}

	else if(key=="Additional information")
	{
	    header.additional_information = values[0];
	}

        else
        {
            // Unknown SDD field
            std::cerr << "Unknown SDD header field entry = " << key << "\n";

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


void SDDparser::printHeaderSummary(std::ostream& out) const {
// MODIFY TO PRINT HEADER FIELDS
	out << "This is a test that SDDparser reads the SDD header ";
}

