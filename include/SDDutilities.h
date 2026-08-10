#ifndef SDD_UTILITIES_H
#define SDD_UTILITIES_H

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cctype>

// Trimming all white spaces in a line of an SDD file
inline std::string trim(const std::string& s)
{
    size_t first = s.find_first_not_of(" \t\r\n");

    if(first == std::string::npos)
        return "";

    size_t last = s.find_last_not_of(" \t\r\n");

    return s.substr(first,last-first+1);
}

// Splitting all SDD file lines by a delimiter of choice
inline std::vector<std::string> split(const std::string& line,char delimiter)
{
    std::vector<std::string> tokens;

    std::stringstream ss(line);

    std::string item;

    while(std::getline(ss,item,delimiter))
        tokens.push_back(trim(item));

    return tokens;
}

// Useful for header field comparison converting field names to lower case
inline std::string toLower(const std::string& input)
{
    std::string result = input;

    std::transform(
        result.begin(),
        result.end(),
        result.begin(),
        [](unsigned char c)
        {
            return static_cast<char>(std::tolower(c));
        }
    );

    return result;
}

// Removes all whitespaces in header field names and changes them to lower case for comparison.
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
	    continue;

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
	    continue;

        output.push_back(std::stoi(v));
    }


    return output;
}

// --------------------------------------------------------------------- //
// Helper functions to interpret codes in SDD header fields
// --------------------------------------------------------------------- //

// Converting integer specified in "Source type" header field to its meaning.
inline std::string sourceTypeMeaning(int sourceType)
{
    switch(sourceType)
    {
        case 1:
            return "Specified one or multiple monoenergetic particles as incident particle source\n";

        case 2:
            return "Specified one or multiple particles with energy distributions as incident particle source -- Ensure 'Energy distribution' is specified\n";

        case 3:
            return "Specified a phase space incident source\n";

        case 4:
            return "Specified a galactic cosmic ray incident source\n";

        default:
            return "Specified an unknown or no incident source type\n";
    }
}

// Converting only certain common incident particle codes, otherwise
// See Particle Data Group website for particle codes. 
// Care about electrons(11), positrons(-11), photons (22), protons (2212),
// neutrons (2112) and alpha particles (1000020040).
inline std::string incidentParticlesMeaning(const std::vector<int>& incidentParticles)
{

// Check if no incident particles specified
    if (incidentParticles.empty())
    {
        return "No incident particles specified.\n";
    }
    
    std::string result = "";

    for (const auto& id : incidentParticles)
    {
    	switch(id)	
        {
       	    case 11:
	    	result += "Specified ELECTRONS as the incident particles\n";
		break;

	    case -11:
    	        result += "Specified POSITRONS as the incident particles\n";
		break;

	    case 22:
	        result += "Specified PHOTONS as the incident particles\n";
		break;

	    case 2212:
	        result += "Specified PROTONS as the incident particles\n";
		break;

	    case 2112:
	        result += "Specified NEUTRONS as the incident particles\n";
		break;

	    case 1000020040:
	        result += "Specified ALPHA PARTICLES as the incident particles\n";
		break;
	    default:
	        result += "See Particle Data Group Particle ID specifications\n";
        	break;
	}
    }

	return result;
}

// Function to report if a single-track exposure, dose or fluence was specified and at what Gy or particles/um^2
inline std::string doseOrFluenceMeaning(const std::vector<double>& doseOrFluenceVec)
{
    if (doseOrFluenceVec.empty())
    {
 	return "No dose or fluence entries specified\n";
    }

    std::string result = "";

// Check if first value specified in the field is 0, 1, or 2 then handle each case.
    if (doseOrFluenceVec[0] == static_cast<double>(0))
    {
        result += "Single-track irradiation specified -- no dose or fluence specifications\n";
    }
    else if(doseOrFluenceVec[0] == static_cast<double>(1))
    {
        result += "Specified a delivered dose of " + std::to_string(doseOrFluenceVec[1]) + " Gy was specified\n";
    }
    else if(doseOrFluenceVec[0] == static_cast<double>(2))
    {
	result += "Specified a fluence of " + std::to_string(doseOrFluenceVec[1]) + "particles per square micrometer\n";    
    }
// If first value is not 0, 1, or 2 return unknown value specified
    else
    {
	result += "Unknown value specified in 'Dose or fluence' field, the first entry should be either '0', '1', or '2'\n";
    }

    return result;
} 


// Interpreting numerical field in the Volume header field.
inline std::string volumesMeaning(const std::vector<double>& volumesVec)
{

    std::string result = "";
// Volumes field can have either 7 or 14 data entries, 7 if the cell and nuclear
// have the same dimensions, 14 to specify the cell dimensions (1-7) and
// the nucleus dimensions (8-14)
    if (volumesVec.size() == 7)  // cell and nucleus have same specifications
    {
	switch(static_cast<int>(volumesVec[0]))
	{
	    case 0: // 0 = cell/nucleus in shape of box
	        result += "The cell and cell nucleus are modeled as BOXES with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers\n";
	        result += "The corresponding Euler rotations of the cell and cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ")\n";
		break;
	    case 1: // 1 = cell/nucleus in shape of ellipsoid
		result += "The cell and cell nucleus are modeled as ELLIPSOIDS with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers\n";
		result += "The corresponding Euler rotations of the cell and cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ")\n";
		break;
	    case 2: // 2 = cell/nucleus in shape of cylinder
		result += "The cell and cell nucleus are modeled as CYLINDERS with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers\n";
		result += "The corresponding Euler rotations of the cell and cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ")\n";
		break;
	    default:
		result += "Invalid 'Volumes' first field entry, please specify either '0' (box), '1' (ellipsoid), '2' (cylinder)";
		break;
	}

    }
// First 7 entries are for the cell, next 7 entries are for the nucleus/scoring volume.
    else if (volumesVec.size() == 14)
    {

	switch(static_cast<int>(volumesVec[0])) //First 7 entries for the cell 
	{
	    case 0: // 0 = cell in shape of box
	        result += "The cell is modeled as a BOX with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers\n";
	        result += "The corresponding Euler rotations of the cell are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ")\n";
		break;
	    case 1: // 1 = cell in shape of ellipsoid
		result += "The cell is modeled as an ELLIPSOID with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers\n";
		result += "The corresponding Euler rotations of the cell are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ")\n";
		break;
	    case 2: // 2 = cell in shape of cylinder
		result += "The cell is modeled as a CYLINDER with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[1]) + ", " + std::to_string(2*volumesVec[2]) + ", " + std::to_string(2*volumesVec[3]) + ") micrometers\n";
		result += "The corresponding Euler rotations of the cell are (phi, theta, psi) = (" + std::to_string(volumesVec[4]) + ", " + std::to_string(volumesVec[5]) + ", " + std::to_string(volumesVec[6]) + ")\n";
		break;
	    default:
		result += "Invalid 'Volumes' first field entry, please specify either '0' (box), '1' (ellipsoid), '2' (cylinder)";
		break;
	}

// Second set of switch cases for the nucleus in fields 8 to 14
	switch(static_cast<int>(volumesVec[7]))
	{
	    case 0: // 0 = nucleus in shape of box
	        result += "The cell nucleus is modeled as a BOX with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[8]) + ", " + std::to_string(2*volumesVec[9]) + ", " + std::to_string(2*volumesVec[10]) + ") micrometers\n";
	        result += "The corresponding Euler rotations of the cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[11]) + ", " + std::to_string(volumesVec[12]) + ", " + std::to_string(volumesVec[13]) + ")\n";
		break;
	    case 1: // 1 = cell/nucleus in shape of ellipsoid
		result += "The cell nucleus is modeled as an ELLIPSOID with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[8]) + ", " + std::to_string(2*volumesVec[9]) + ", " + std::to_string(2*volumesVec[10]) + ") micrometers\n";
		result += "The corresponding Euler rotations of the cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[11]) + ", " + std::to_string(volumesVec[12]) + ", " + std::to_string(volumesVec[13]) + ")\n";
		break;
	    case 2: // 2 = cell/nucleus in shape of cylinder
		result += "The cell nucleus is modeled as a CYLINDER with side lengths (X, Y, Z) = (" + std::to_string(2*volumesVec[8]) + ", " + std::to_string(2*volumesVec[9]) + ", " + std::to_string(2*volumesVec[10]) + ") micrometers\n";
		result += "The corresponding Euler rotations of the cell nucleus are (phi, theta, psi) = (" + std::to_string(volumesVec[11]) + ", " + std::to_string(volumesVec[12]) + ", " + std::to_string(volumesVec[13]) + ")\n";
		break;
	    default:
		result += "Invalid 'Volumes' eighth field entry, please specify either '0' (box), '1' (ellipsoid), '2' (cylinder)";
		break;
	}

    }

    else
    {
	result += "Incorrect number of fields specified for the Volumes section (please ensure either 7 or 14 fields are specified";
    }

    return result;
}


inline std::string cellCyclePhaseMeaning(const std::vector<double>& cellCycleVec)
{
    if (cellCycleVec.empty())
    {
	return "No Cell cycle header entry specified";
    }

    std::string result = "";

    switch(static_cast<int>(cellCycleVec[0]))
    {
	case 1: 
	    result += "Cells are " + std::to_string(100*cellCycleVec[1]) + "% through the G0 phase";
		break;
	case 2: 
	    result += "Cells are " + std::to_string(100*cellCycleVec[1]) + "% through the G1 phase";
		break;
	case 3: 
	    result += "Cells are " + std::to_string(100*cellCycleVec[1]) + "% through the S phase";
		break;
	case 4: 
	    result += "Cells are " + std::to_string(100*cellCycleVec[1]) + "% through the G2 phase";
		break;
	case 5: 
	    result += "Cells are " + std::to_string(100*cellCycleVec[1]) + "% through the M phase";
		break;
	default:
	    result += "Cell cycle phase unknown, if desired, specify: 1 (G0), 2 (G1), 3 (S), 4 (G2), or 5 (M).\n";
		break;
    }

	return result;
}

inline std::string dnaStructureMeaning(const std::vector<int>& dnaStructureVec)
{
    if (dnaStructureVec.empty())
    {
	return "No DNA structure information specified.\n";
    }

    std::string result = "";

    switch(dnaStructureVec[0])
    {
	case 0:
	    result += "DNA structured in a 'Whole Nucleus' arrangement\n";
	    break;
    	case 1:
	    result += "DNA structured as a 'Heterochromatin' region\n";
	    break;
        case 2:
	    result += "DNA structured as a 'Euchromatin' region\n";
	    break;
        case 3:
	    result += "DNA structured in a 'Mixed Heterochromatin and Euchromatin' region\n";
	    break;
        case 4:
	    result += "DNA structured as a 'Single DNA fiber'\n";
	    break;
        case 5:
	    result += "DNA structured as DNA wrapped around a single histone\n";
	    break;
        case 6:
	    result += "DNA structured as a 'DNA plasmid'\n";
	    break;
        case 7:
	    result += "DNA structured as a 'Simple circular section'\n";
	    break;
	case 8:
	    result += "DNA structured as a 'Simple straight section'\n";
	    break;
        default:
 	    result += "No valid DNA structure entry provided.";
	    break;
    }

	return result;
}


inline std::string proliferationStatusMeaning(int status)
{
    if (status == 0)
    {
	return "Cell(s) are in a quiescent state";
    }
    else if (status == 1)
    {
	return "Cell(s) are in a proliferating state";
    }
    else
    {
   	return "Proliferation Status entry invalid";
    }
}


inline std::string damageDefinitionMeaning(const std::vector<double>& damageDefinitionVec)
{
    if (damageDefinitionVec.empty())
    {
	return "No Damage definition information specified";
    }

    std::string result = "";

    switch(static_cast<int>(damageDefinitionVec[0]))
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

    if (static_cast<int>(damageDefinitionVec[1]) == 0)
    {
	result += "Backbone lesions that are considered as DSBs are within a distance of " + std::to_string(damageDefinitionVec[2]) + " bp.\n";

	if (static_cast<int>(damageDefinitionVec[3]) == -1)
	{	
	    result += "Base lesions are not being scored.\n";
    	}

    	else if (static_cast<int>(damageDefinitionVec[3]) >= 0)
    	{
            result += "Base damages a distance of " + std::to_string(damageDefinitionVec[3]) + " bp beyond the outer backbone damages are stored in the same site.\n";
    	}

    	else 
    	{
            result += "Invalid Damage definition fielf entry encountered in field 4, please specify either a positive integer, 0, or -1.\n";
    	}
    }

    else if (static_cast<int>(damageDefinitionVec[1]) == 1)
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
            result += "Invalid Damage definition fielf entry encountered in field 4, please specify either a positive integer, 0, or -1.\n";
    	}	
    }

    else
    {
	result += "Invalid second field entry for Damage definition, please choose either 0 (bp) or 1 (nm).\n";
    }


    result += "The lower energy threshold specified to induce a strand break or base damage is: " + std::to_string(damageDefinitionVec[4]) + " eV.\n";

    return result;
}







#endif
