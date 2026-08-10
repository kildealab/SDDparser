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

// FIX THIS
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


void SDDparser::printHeaderSummary(std::ostream& output) const {

// First section is the radiation information
	output << "--------------------------- Incident Radiation Information "
	       << "---------------------------\n\n";

// Source summary - prints the string entry
	output << "Radiation source: " << header.source << "\n";
	
// Source Type summary
	output << "Source type: " << header.source_type << " ("
	<< sourceTypeMeaning(header.source_type) << ")\n";

// Incident Particles summary
	output << incidentParticlesMeaning(header.incident_particles) << "\n" <<
	"Each specified incident particle had the following fluence fractions in that order: ";
	for (size_t i = 0; i < header.particle_fraction.size(); i++)
	{
	    output << header.particle_fraction[i];
	}
	output << "\n";

// Dose or fluence summary
	output << doseOrFluenceMeaning(header.dose_or_fluence) << "\n"; 


// Second section is the target information
	output << "-------------------------- Radiation Target Information" <<
		  "--------------------------\n";
// Irradiation target
	output << "Radiation incident on: " << header.irradiation_target << "\n";

// Summary of target shapes and volumes
	output << volumesMeaning(header.volumes) << "\n" ;

// Summary of chromosome number and sizes, DNA density, cell cycle phase, proliferation status.
	output << "The number of chromosomes specified were: " << header.chromosome_sizes[0] << ". Subsequent field entries are the chromosome sizes in units of Mega base pairs (Mbp),\nwith an average DNA density of " << header.DNA_density << " Mbp per cubic micrometer.\n";
	output << cellCyclePhaseMeaning(header.cell_cycle_phase) << "\n" ;
	output << "The DNA is structured as: " << header.DNA_structure[0] << " (" << dnaStructureMeaning(header.DNA_structure) << ")\n";
	output << "The cell Proliferation status is: " << header.proliferation_status << " (" << proliferationStatusMeaning(header.proliferation_status) << ") \n\n";

//Third section is the DNA damage and exposure information
	output << "------------------------- DNA Damage Information" <<
		  "-------------------------\n";

// Summary of the Damage definition and damage and primary count header fields
	output << damageDefinitionMeaning(header.damage_definition) << "\n";
	output << "The number of distinct damage lesions scored is " << header.damage_and_primary_count[0] << ", as a result of " << header.damage_and_primary_count[1] << " primary particles simulated.\n";

}
void SDDparser::printExposureSummary(std::ostream& output) const {
	
	output << "Test for printExposureSummary\n";
}

void SDDparser::printSummary(std::ostream& output) const
{
    output << "========================================\n";
    output << "             SDD FILE SUMMARY\n";
    output << "========================================\n\n";

    printHeaderSummary(output);

    output << "\n";

    printExposureSummary(output);

    output << "\n";
    output << "========================================\n";
}
