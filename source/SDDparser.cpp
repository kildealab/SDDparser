#include "SDDparser.h"
#include "SDDutilities.h"

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

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

    if(!parseDamageEntries(file))
    {
    	std::cerr << "ERROR: SDD data entry parsing failed.\n";
    	return false;
    }

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

bool SDDparser::parseDamageEntries(std::ifstream& file)
{

    std::string line;

    std::size_t expectedFields = 0;

    int currentExposureID = 0;

    for(int enabled : header.data_entries)
    {
	if(enabled == 1)
        expectedFields++;
    }

    while(std::getline(file,line))
    {
	line = trim(line);

        if(line.empty())
            continue;

        DamageEntry damage;

        damage.rawLine = line;

	//--------------------------------------------------
        // Split the data entry into Fields 1-7
        //--------------------------------------------------

        std::vector<std::string> fields = split(line, ';');

        if(fields.size() != expectedFields)
        {
            std::cerr << "ERROR: SDD data entry contains "
            << fields.size()
            << " fields, but the header specifies "
            << expectedFields
            << " fields.\n";

    	    return false;
	}

	
	std::size_t fieldIndex = 0;
 	//--------------------------------------------------
        // Field 1: Classification
        //--------------------------------------------------

        if(header.data_entries.size() > 0 && header.data_entries[0] == 1)
        {

 	    if(fieldIndex >= fields.size())
    	    {
        	std::cerr << "ERROR: Missing SDD data entry Field 1.\n";
        	return false;
    	    }

	    std::vector<std::string> values = split(fields[fieldIndex], ',');

            if(values.size() >= 1 && !trim(values[0]).empty())
            {

	        int exposureMarker = std::stoi(trim(values[0]));
		damage.classification.exposureMarker = exposureMarker;

// A value of 2 indicates the beginning of a new exposure. If SDD data field 1
// does not start with a 2, force it to be considered a new exposure anyway.      	
        	if(exposureMarker == 2 || currentExposureID == 0)
        	{
            	    currentExposureID++;

            	    Exposure newExposure;
            	    newExposure.exposureID = currentExposureID;

            	    exposures.push_back(newExposure);
        	}
	    }

            if(values.size() >= 2 && !trim(values[1]).empty())
            {
	        damage.classification.eventID = std::stoi(trim(values[1]));
	    }

	    fieldIndex++;
        }


	//--------------------------------------------------
        // Field 2: Spatial Position
        //--------------------------------------------------

        if(header.data_entries.size() > 1 && header.data_entries[1] == 1)
        {

 	    if(fieldIndex >= fields.size())
    	    {
        	std::cerr << "ERROR: Missing SDD data entry Field 2.\n";
        	return false;
    	    }
	    std::vector<std::string> fieldEntries;

	    std::vector<std::string> subfields = split(fields[fieldIndex], '/');
	    if (subfields.size() == 3)
	    {
		for (const std::string& subfield: subfields)
		{
       	            std::vector<std::string> values = split(subfield, ',');
		    fieldEntries.insert(fieldEntries.end(), values.begin(), values.end());
	        }

		if (fieldEntries.size() != 9)
		{
		    std::cerr << "ERROR: Incorrect number of entries for Data field 2.\n";
    		    return false;
		}

		else
		{
		damage.position.x = std::stod(trim(fieldEntries[0]));
		damage.position.y = std::stod(trim(fieldEntries[1]));
		damage.position.z = std::stod(trim(fieldEntries[2]));
         	damage.position.x_max = std::stod(trim(fieldEntries[3]));
		damage.position.y_max = std::stod(trim(fieldEntries[4]));
		damage.position.z_max = std::stod(trim(fieldEntries[5]));
         	damage.position.x_min = std::stod(trim(fieldEntries[6]));
		damage.position.y_min = std::stod(trim(fieldEntries[7]));
		damage.position.z_min = std::stod(trim(fieldEntries[8]));
         	}

	    }

	    else if (subfields.size() == 1)
	    {
		std::vector<std::string> values = split(fields[fieldIndex], ',');
		damage.position.x = std::stod(trim(values[0]));
		damage.position.y = std::stod(trim(values[1]));
		damage.position.z = std::stod(trim(values[2]));
            }
	    else
	    {
		std::cerr << "ERROR: Incorrect number of entries for Data field 2.\n";
		return false;
	    }
	    
    	    fieldIndex++;
        }

	//--------------------------------------------------
        // Field 3: Chromosome IDs
        //--------------------------------------------------

        if(header.data_entries.size() > 2 && header.data_entries[2] == 1)
        {

 	    if(fieldIndex >= fields.size())
    	    {
        	std::cerr << "ERROR: Missing SDD data entry Field 3.\n";
        	return false;
    	    }

            std::vector<std::string> values = split(fields[fieldIndex], ',');

            if(values.size() >= 1 && !trim(values[0]).empty())
            {    
		damage.chromosomeID.dnaStructure = std::stoi(trim(values[0]));
	    }
            if(values.size() >= 2 && !trim(values[1]).empty())
            {
		damage.chromosomeID.chromosomeNumber = std::stoi(trim(values[1]));
	    }
            if(values.size() >= 3 && !trim(values[2]).empty())
            {
		damage.chromosomeID.chromatidNumber = std::stoi(trim(values[2]));
	    }
            if(values.size() >= 4 && !trim(values[3]).empty())
            {
	        damage.chromosomeID.chromosomeArm = std::stoi(trim(values[3]));
            }	

	    fieldIndex++;
	}

	//--------------------------------------------------
        // Field 4: Chromosome Position
        //--------------------------------------------------

        if(header.data_entries.size() > 3 && header.data_entries[3] == 1)
        {

 	    if(fieldIndex >= fields.size())
    	    {
        	std::cerr << "ERROR: Missing SDD data entry Field 4.\n";
        	return false;
    	    }

            double value = std::stod(trim(fields[fieldIndex]));

            damage.chromosomePosition.position = value;

            if(value < 1.0)
                damage.chromosomePosition.isFractional = true;
            else
                damage.chromosomePosition.isFractional = false;

	    fieldIndex++;
        }



 	//--------------------------------------------------
        // Field 5: Cause
        //--------------------------------------------------

        if(header.data_entries.size() > 4 && header.data_entries[4] == 1)
        {

 	    if(fieldIndex >= fields.size())
    	    {
        	std::cerr << "ERROR: Missing SDD data entry Field 5.\n";
        	return false;
    	    }

            std::vector<std::string> values = split(fields[fieldIndex], ',');

            if(values.size() >= 1 && !values[0].empty())
                damage.damageCause.cause =
                    std::stoi(trim(values[0]));

            if(values.size() >= 2 && !values[1].empty())
                damage.damageCause.numDirectDamages =
                    std::stoi(trim(values[1]));

            if(values.size() >= 3 && !values[2].empty())
                damage.damageCause.numIndirectDamages =
                    std::stoi(trim(values[2]));
        
	    fieldIndex++;
	}


        //--------------------------------------------------
        // Field 6: Damage Types
        //--------------------------------------------------

        if(header.data_entries.size() > 5 && header.data_entries[5] == 1)
        {

 	    if(fieldIndex >= fields.size())
    	    {
        	std::cerr << "ERROR: Missing SDD data entry Field 6.\n";
        	return false;
    	    }

            std::vector<std::string> values = split(fields[fieldIndex], ',');

            if(values.size() >= 1 && !values[0].empty())
                damage.damageType.numBaseDamages =
                    std::stoi(trim(values[0]));

            if(values.size() >= 2 && !values[1].empty())
                damage.damageType.numSingleBackboneBreaks =
                    std::stoi(trim(values[1]));

            if(values.size() >= 3 && !values[2].empty())
                damage.damageType.presenceOfDSB =
                    std::stoi(trim(values[2]));
        
	    fieldIndex++;
	}


        //--------------------------------------------------
        // Field 7: Full Break Specification
        //--------------------------------------------------

        if(header.data_entries.size() > 6 && header.data_entries[6] == 1)
        {

 	    if(fieldIndex >= fields.size())
    	    {
        	std::cerr << "ERROR: Missing SDD data entry Field 7.\n";
        	return false;
    	    }

            std::vector<std::string> groups = split(fields[fieldIndex], '/');

            for(const auto& group : groups)
            {
                if(trim(group).empty())
                {
                    continue;
                }
		
		// Each group consists of: STRAND, BASE, baseDamageType
        	std::vector<std::string> values = split(group, ',');
    		
		if(values.size() != 3)
        	{
            	    std::cerr << "ERROR: Invalid Field 7 group: [" << group << "]\n";
            	    return false;
        	}
         

	        int strand = std::stoi(trim(values[0]));

        	int base = std::stoi(trim(values[1]));

        	int baseDamageType = std::stoi(trim(values[2]));

		// Validate STRAND
        	if(strand < 1 || strand > 4)
        	{
            	    std::cerr << "ERROR: Invalid Field 7 strand value: " << strand << "\n";

            	    return false;
        	}

        	// Validate base damage type
        	if(baseDamageType < 0 || baseDamageType > 3)
        	{
            	    std::cerr << "ERROR: Invalid Field 7 damage type: " << baseDamageType << "\n";

            	    return false;
        	}


		damage.fullBreakSpec.strand.push_back(strand);
        	damage.fullBreakSpec.base.push_back(base);
        	damage.fullBreakSpec.baseDamageType.push_back(baseDamageType);

	   }
        
	    fieldIndex++;
	}


        //--------------------------------------------------
        // Add damage to the appropriate exposure
        //--------------------------------------------------

        
	if(currentExposureID <= 0 || currentExposureID > static_cast<int>(exposures.size()))
	{
    	    std::cerr << "ERROR: Invalid current exposure ID.\n";
   	    return false;
	}

	exposures[currentExposureID - 1].damages.push_back(damage);

    }

    return true;
}


const Header& SDDparser::getHeader() const
{
    return header;
}


const std::vector<Exposure>& SDDparser::getExposures() const 
{
    return exposures;
}


std::map<int, ChromosomeDamageSummary> SDDparser::summarizeChromosomeDamage(const Exposure& exposure) const
{
    std::map<int, ChromosomeDamageSummary> summary;

    for(const auto& damage : exposure.damages)
    {
        const auto& chromosome = damage.chromosomeID;
        const auto& damageType = damage.damageType;

        // Only strands that have an associated chromosome number 
	// (non-zero integer) will be valid.
        if(chromosome.chromosomeNumber <= 0)
        {
            continue;
        }

        int chromosomeNumber = chromosome.chromosomeNumber;

        ChromosomeDamageSummary& chromosomeSummary = summary[chromosomeNumber];

        chromosomeSummary.dnaStructure = chromosome.dnaStructure;

        chromosomeSummary.chromosomeNumber = chromosomeNumber;

        chromosomeSummary.numBaseDamages += damageType.numBaseDamages;

        chromosomeSummary.numSingleStrandBreaks += damageType.numSingleBackboneBreaks;

        chromosomeSummary.numDSBs += damageType.presenceOfDSB;
    }

    return summary;
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

// Summary of microenvironment and time for which the chemistry simulation ends.
	output << microenvironmentMeaning(header.microenvironment) << "\n";
	output << timeMeaning(header.time) << "\n";

// Summary of Data field entries
	output << dataEntriesMeaning(header.data_entries) << "\n";

}


void SDDparser::printExposureSummary(std::ostream& output) const {
	output << "-----------------------------------------------/n";	
	output << "------------------------ Chromosome Damages" <<
		  "------------------------\n";
	output << "-----------------------------------------------/n";

	output << "\nNumber of exposures: " << exposures.size() << "\n";


    for(const auto& exposure : exposures)
    {
        output << "\n----------------------------------------\n";
        output << "Exposure " << exposure.exposureID << "\n";
        output << "----------------------------------------\n";

        output << "Number of damage entries: "
            << exposure.damages.size()
            << "\n";

        auto chromosomeSummary =
            summarizeChromosomeDamage(exposure);

        if(chromosomeSummary.empty())
        {
            output << "\nNo chromosome-associated damage entries found.\n";
            continue;
        }

        output << "\nDamage by chromosome:\n\n";

	int totalBaseDamages = 0;
	int totalSingleStrandBreaks = 0;
	int totalDSBs = 0;

        for(const auto& [chromosomeNumber, summary] :
            chromosomeSummary)
        {
	    totalBaseDamages += summary.numBaseDamages;
	    totalSingleStrandBreaks += summary.numSingleStrandBreaks;
	    totalDSBs += summary.numDSBs;

            output << "Chromosome " << chromosomeNumber << "\n";

            output << "  Base damages: " << summary.numBaseDamages << "\n";

            output << "  Single-strand breaks: " << summary.numSingleStrandBreaks << "\n";

            output << "  Double-strand breaks: " << summary.numDSBs << "\n";

	}

	output << "\nExposure totals:\n";

	output << "  Base damages: " << totalBaseDamages << "\n";

	output << "  Single-strand breaks: " << totalSingleStrandBreaks << "\n";

	output << "  Double-strand breaks: " << totalDSBs << "\n";

    }
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
