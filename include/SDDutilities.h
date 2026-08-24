#ifndef SDD_UTILITIES_H
#define SDD_UTILITIES_H

#include <cstddef>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cctype>

// -------------------------------------------------- //
// ----------- GENERIC HELPER FUNCTIONS ------------- //
// -------------------------------------------------- //

// Trimming all white spaces in a line of an SDD file
inline std::string trim(const std::string& s)
{
    size_t first = s.find_first_not_of(" \t\r\n");

    if(first == std::string::npos)
    {
        return "";
    }

    size_t last = s.find_last_not_of(" \t\r\n");

    return s.substr(first,last-first+1);
}

// Splitting all SDD file lines by a delimiter of choice ',', ';', '/'.
inline std::vector<std::string> split(const std::string& line,char delimiter)
{
    std::vector<std::string> tokens;

    std::stringstream ss(line);

    std::string item;

    while(std::getline(ss,item,delimiter))
    {
	tokens.push_back(trim(item));
    }

    return tokens;
}

// Useful for header field comparison converting field names to lower case
inline std::string toLower(const std::string& input)
{
    std::string result = input;

    std::transform(result.begin(), 
	    result.end(), 
	    result.begin(), 
	    [](unsigned char c)
        {
            return static_cast<char>(std::tolower(c));
        }
    );
    return result;
}

// Removes all whitespaces in header field names and changes names to lower case for comparison.
inline std::string normalizeHeaderKey(const std::string& input)
{
    return toLower(trim(input));
}

// Helper function to parse vectors consisting of doubles
inline std::vector<double> parseDoubleList(const std::vector<std::string>& values)
{
    std::vector<double> output;

    for(auto& v : values)
    {
        if (v.empty())
	{    
	     continue;
	}
	output.push_back(std::stod(v));
    }

    return output;
}

// Helper function to parse vectors consiting of ints
inline std::vector<int> parseIntList(const std::vector<std::string>& values)
{
    std::vector<int> output;

    for(auto& v : values)
    {
	if (v.empty())
	{
	    continue;
	}
        output.push_back(std::stoi(v));
    }

    return output;
}

// --------------------------------------------------------------------- //
// ------ HELPER FUNCTIONS TO INTERPRET SDD HEADER FIELD ENTRIES ------- //
// --------------------------------------------------------------------- //

// Return meaning of numeric entry in 'Source type' header field
inline std::string sourceTypeMeaning(int sourceType)
{
    switch(sourceType)
    {
        case 1:
            return "Specified one or multiple monoenergetic particles as incident particle source";

        case 2:
            return "Specified one or multiple particles with energy distributions as incident particle source -- Ensure 'Energy distribution' is specified";

        case 3:
            return "Specified a phase space incident source";

        case 4:
            return "Specified a galactic cosmic ray incident source";

        default:
            return "Specified an unknown or no incident source type";
    }
}

// Converting only certain common incident particle codes, otherwise
// See Particle Data Group website for particle codes. 
// Most common are: electrons(11), positrons(-11), photons (22), 
// protons (2212), neutrons (2112) and alpha particles (1000020040).
// Function to be called in a for loop through each incident particle listed 
// in the header to be summarized with the particles associated particle fluence
// fraction.
inline std::string incidentParticlesMeaning(int incidentParticle)
{

    std::string result = "";

    switch(incidentParticle)	
    {
       	case 11:
   	    result += "ELECTRONS";
	    break;

        case -11:
            result += "POSITRONS";
       	    break;

	case 22:
	    result += "PHOTONS";
	    break;

	case 2212:
	    result += "PROTONS";
	    break;

	case 2112:
	    result += "NEUTRONS";
	    break;

	case 1000020040:
	    result += "ALPHA PARTICLES";
	    break;

	default:
	    result += "See Particle Data Group Particle ID specifications\n";
            break;
    }
	return result;
}

// Function to report if a single-track exposure, dose or fluence was specified and at what Gy or particles/um^2
inline std::string doseOrFluenceMeaning(const std::vector<double>& doseOrFluenceVec)
{
    if (doseOrFluenceVec.empty())
    {
 	return "None specified.";
    }

    std::string result = "";

    if (doseOrFluenceVec[0] == static_cast<double>(0)) 		// If 0 is specified = a single-track irradiation
    {
        result += "Single-track irradiation.";
    }
    else if(doseOrFluenceVec[0] == static_cast<double>(1))	// If 1 is specified = a delivered dose
    {
        result += "Dose = " + std::to_string(doseOrFluenceVec[1]) + " Gy.";
    }
    else if(doseOrFluenceVec[0] == static_cast<double>(2))	// If 2 is specified = a fluence
    {
	result += "Fluence = " + std::to_string(doseOrFluenceVec[1]) + "particles/um^2.";    
    }
    else 							// If first value is not 0, 1, or 2 return unknown value specified
    {
	result += "Unknown value specified.";
    }
    return result;
} 


// Interpreting numerical field in the Volume header field.
inline std::string volumesMeaning(const std::vector<double>& volumesVec)
{
    std::string result = "";
// Volumes field can have either 7 or 14 data entries, 7 if the cell and nucleus
// have the same dimensions, 14 to specify the cell dimensions (1-7) and
// the nucleus dimensions (8-14).
    if (volumesVec.size() == 7)  // cell and nucleus have same specifications
    {
	switch(static_cast<int>(volumesVec[0]))
	{
	    case 0: // 0 = cell/nucleus in shape of box, return box dimensions
	        result += "The cell and cell nucleus are modeled as BOXES with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers.\n";
	        result += "The corresponding Euler rotations of the cell and cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ").\n";
		break;
	    case 1: // 1 = cell/nucleus in shape of ellipsoid, return ellipsoid dimensions
		result += "The cell and cell nucleus are modeled as ELLIPSOIDS with axes lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers.\n";
		result += "The corresponding Euler rotations of the cell and cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ").\n";
		break;
	    case 2: // 2 = cell/nucleus in shape of cylinder, return cylinder dimensions
		result += "The cell and cell nucleus are modeled as CYLINDERS with X and Y axes lengths and Z height (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers.\n";
		result += "The corresponding Euler rotations of the cell and cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ").\n";
		break;
	    default: // If not 0, 1, or 2, return error.
		result += "Invalid 'Volumes' first field entry, please specify either '0' (box), '1' (ellipsoid), '2' (cylinder).";
		break;
	}

    }
// First 7 entries are for the cell, next 7 entries are for the nucleus/scoring volume.
    else if (volumesVec.size() == 14)
    {

	switch(static_cast<int>(volumesVec[0])) //First 7 entries for the cell 
	{
	    case 0: // 0 = cell in shape of box
	        result += "The cell is modeled as a BOX with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers.\n";
	        result += "The corresponding Euler rotations of the cell are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ").\n";
		break;
	    case 1: // 1 = cell in shape of ellipsoid
		result += "The cell is modeled as an ELLIPSOID with axes lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers.\n";
		result += "The corresponding Euler rotations of the cell are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ").\n";
		break;
	    case 2: // 2 = cell in shape of cylinder
		result += "The cell is modeled as a CYLINDER with X and Y axes lengths and Z height (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers.\n";
		result += "The corresponding Euler rotations of the cell are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ").\n";
		break;
	    default:
		result += "Invalid 'Volumes' first field entry, please specify either '0' (box), '1' (ellipsoid), '2' (cylinder).";
		break;
	}

// Second set of switch cases for the nucleus in fields 8 to 14
	switch(static_cast<int>(volumesVec[7]))
	{
	    case 0: // 0 = nucleus in shape of box
	        result += "The cell nucleus is modeled as a BOX with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[8]) + ", " + std::to_string(2*volumesVec[9]) + ", " + std::to_string(2*volumesVec[10]) + ") micrometers.\n";
	        result += "The corresponding Euler rotations of the cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[11]) + ", " + std::to_string(volumesVec[12]) + ", " + std::to_string(volumesVec[13]) + ").\n";
		break;
	    case 1: // 1 = cell/nucleus in shape of ellipsoid
		result += "The cell nucleus is modeled as an ELLIPSOID with axes lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[8]) + ", " + std::to_string(2*volumesVec[9]) + ", " + std::to_string(2*volumesVec[10]) + ") micrometers.\n";
		result += "The corresponding Euler rotations of the cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[11]) + ", " + std::to_string(volumesVec[12]) + ", " + std::to_string(volumesVec[13]) + ").\n";
		break;
	    case 2: // 2 = cell/nucleus in shape of cylinder
		result += "The cell nucleus is modeled as a CYLINDER with X and Y axes lengths and Z height (X, Y, Z) = (" + std::to_string(2*volumesVec[8]) + ", " + std::to_string(2*volumesVec[9]) + ", " + std::to_string(2*volumesVec[10]) + ") micrometers.\n";
		result += "The corresponding Euler rotations of the cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[11]) + ", " + std::to_string(volumesVec[12]) + ", " + std::to_string(volumesVec[13]) + ").\n";
		break;
	    default:
		result += "Invalid 'Volumes' eighth field entry, please specify either '0' (box), '1' (ellipsoid), '2' (cylinder).";
		break;
	}

    }

    else 		// If not 7 or 14 entries, exit with error.
    {
	result += "Incorrect number of fields specified for the Volumes section (please ensure either 7 or 14 fields are specified.";
    }

    return result;
}


// Function to return which phase of the cell cycle the cell model is in.
inline std::string cellCyclePhaseMeaning(const std::vector<double>& cellCycleVec)
{
    if (cellCycleVec.empty())			// If no value specified, return no data present.
    {
	return "None specified.";
    }

    std::string result = "";

    switch(static_cast<int>(cellCycleVec[0])) // Look at first integer entry for cell cycle phase
    {
	case 1: // G0 phase
	    result += "G0.";
		break;
	case 2: // G1 phase
	    result += "G1.";
		break;
	case 3: // S phase
	    result += "S.";
		break;
	case 4: // G2 phase
	    result += "G2.";
		break;
	case 5: // M phase
	    result += "M.";
		break;
	default: // if other number specified, exit with warning.
	    result += "Unknown.";
		break;
    }
	return result;
}

// Function to return intepretation of DNA structure entry (value of 0-8)
inline std::string dnaStructureMeaning(const std::vector<int>& dnaStructureVec)
{
    if (dnaStructureVec.empty())		// If none specified, raise warning.
    {
	return "No DNA structure information specified.\n";
    }

    std::string result = "";

    switch(dnaStructureVec[0])
    {
	case 0:
	    result += "DNA structured in a 'Whole Nucleus' arrangement";
	    break;
    	case 1:
	    result += "DNA structured as a 'Heterochromatin' region";
	    break;
        case 2:
	    result += "DNA structured as a 'Euchromatin' region";
	    break;
        case 3:
	    result += "DNA structured in a 'Mixed Heterochromatin and Euchromatin' region";
	    break;
        case 4:
	    result += "DNA structured as a 'Single DNA fiber'";
	    break;
        case 5:
	    result += "DNA structured as DNA wrapped around a single histone";
	    break;
        case 6:
	    result += "DNA structured as a 'DNA plasmid'";
	    break;
        case 7:
	    result += "DNA structured as a 'Simple circular section'";
	    break;
	case 8:
	    result += "DNA structured as a 'Simple straight section'";
	    break;
        default:
 	    result += "No valid DNA structure entry provided";
	    break;
    }

	return result;
}

// Proliferation status of a cell can be 0 = quiescent or 1 = proliferating.
inline std::string proliferationStatusMeaning(const std::vector<std::string>& status)
{
    if (std::stoi(status[0]) == 0)
    {
	return "Cell(s) are in a quiescent state";
    }
    else if (std::stoi(status[0]) == 1)
    {
	return "Cell(s) are in a proliferating state";
    }
    else 	// If not 0 or 1, error.
    {
   	return "Proliferation Status entry invalid";
    }
}

// Handling each field entry according to the Damage Definition header
inline std::string damageDefinitionMeaning(const std::vector<double>& damageDefinitionVec)
{
    if (damageDefinitionVec.empty()) // Inform user if not specified
    {
	return "No Damage definition information specified";
    }

    if (damageDefinitionVec.size() != 5 && damageDefinitionVec.size() != 6) // This header can have either 5 fields or 6 fields (6th is optional).
    {
	return "Incorrect number of entries in 'Damage Definition' header\n";
    }

    std::string result = "";

    switch(static_cast<int>(damageDefinitionVec[0])) // Check first entry (either 0 for direct, or 1 for including chemistry)
    {
	case 0:
	    result += "Damages recorded as a result of direct effects only.\n";
	    break;
	case 1:
	    result += "Damaged recorded as a result of direct and indirect (chemistry) effects.\n";
	    break;
	default:
	    result += "Damage definition first entry invalid, please choose either 0 (direct) or 1 (direct + indirect).\n";
	    break;
    }

    if (static_cast<int>(damageDefinitionVec[1]) == 0) // Check second entry (0 indicates definition is in units of base pairs)
    {
	result += "Backbone lesions that are considered as DSBs are within a distance of " + std::to_string(damageDefinitionVec[2]) + " bp.\n"; // Append to the string the value in the third entry

	if (static_cast<int>(damageDefinitionVec[3]) == -1) // Check fourth entry
	{	
	    result += "Base lesions are not being scored.\n";
    	}

    	else if (static_cast<int>(damageDefinitionVec[3]) >= 0)
    	{
            result += "Base damages a distance of " + std::to_string(damageDefinitionVec[3]) + " bp beyond the outer backbone damages are stored in the same site.\n";
    	}

    	else
    	{
            result += "Invalid Damage definition field entry encountered in field 4, please specify either a positive integer, 0, or -1.\n";
    	}
    }

    else if (static_cast<int>(damageDefinitionVec[1]) == 1) // Check second entry (1 indicates units are in nm)
    {

	result += "Backbone lesions that are considered as DSBs are within a distance of " + std::to_string(damageDefinitionVec[2]) + " nm.\n";

    	if (static_cast<int>(damageDefinitionVec[3]) == -1)
    	{
            result += "Base lesions are not being scored.\n";
    	}

    	else if (static_cast<int>(damageDefinitionVec[3]) >= 0)
        {
            result += "Base damages a distance of " + std::to_string(damageDefinitionVec[3]) + " nm beyond the outer backbone damages are stored in the same site.\n";
    	}

    	else 
    	{
            result += "Invalid Damage definition field entry encountered in field 4, please specify either a positive integer, 0, or -1.\n";
    	}	
    }

    else
    {
	result += "Invalid second field entry for Damage definition, please choose either 0 (bp) or 1 (nm).\n";
    }

    result += "The lower energy threshold specified to induce a strand break or base damage is: " + std::to_string(damageDefinitionVec[4]) + " eV.\n"; // Append 5th field entry

    return result;
}

// Header field to specify temperature and oxygen conditions
inline std::string microenvironmentMeaning(const std::vector<double>& microenvVec)
{
    if (microenvVec.empty()) // If none, assume default parameters
    {
	return "No microenvironment information submitted. Assumed the temperature is 25 degrees Celsius and normoxic conditions.\n";
    }

    else if (microenvVec.size() == 2) // If 2 fields specified, first is temperature, second is oxygen level
    {
	return "Microenvironment temperature is: " + std::to_string(microenvVec[0]) + " degrees Celsius. The microenvironment oxygen concentration is: " + std::to_string(microenvVec[1]) + " molarity (M).\n";
    }

    else // If incorrect number of entries, indicate.
    {
	return "Incorrect number of 'Microenvironment' field entries specified - Please specify two floats.\n";
    }
}

// Particle simulation time depending on if chemistry simulations are involved
inline std::string timeMeaning(double t)
{
    if (t == 0) // Only direct physical interactions considered
    {
	return "Time of: " + std::to_string(t) + " seconds specified, indicating simulations only consider direct physics interactions, no chemical interactions.\n";
    }

    else if (t > 0) // Includes chemistry interactions.
    {
	return "Time of: " + std::to_string(t) + "seconds from the time when the source particle was simulated to the time at which the chemistry simulation ends.\n";
    }

    else
    {
	return "Invalid 'Time' field entry encountered, please enter a value greater than or equal to 0.\n";
    }
}

// Essential to interpret the SDD data entries, check the order of 0s and 1s, should have 14 header entries. 
inline std::string dataEntriesMeaning(const std::vector<int>& dataEntriesVec)
{
    if (dataEntriesVec.empty())
    {
	return "Data entries required to interpret exposure data. Please refer to the official SDD paper to fill these fields.";
    }

    if (dataEntriesVec.size() != 14)
    {
	return "Incorrect number of entries for 'Data entries' field, please specify '14' fields.\n";
    }

    std::string result = "The following data field entries have been specified: ";

    bool first_entry = true;

    for (size_t i = 0; i < dataEntriesVec.size(); i++)
    {
	if (dataEntriesVec[i] == 1) // If entry is a 1, indicate that field is being processed in the data block
	{
	    if(!first_entry) // Help formatting by removing trailing commas
	    {
	    	result += ", ";
	    }

	    result += "Field " + std::to_string(i+1);

	    first_entry = false;
	}
	else if (dataEntriesVec[i] == 0) // If entry is 0, skip that field in the data block
	{
	    continue;
	}

	else // If not 0 or 1, error
	{
	    return "Incorrect field entry encountered: " + std::to_string(dataEntriesVec[i]) + ". Please specify either with '0' or '1'.\n";
	}
    }
    result += ".\n";
    return result;
}


#endif
