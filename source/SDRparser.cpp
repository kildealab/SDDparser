#include "SDRparser.h"
#include "utilities.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

bool SDRparser::parseFile(const std::string& filename)
{
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error: Could not open SDR file: "
                  << filename << std::endl;

        return false;
    }

    // Clear previous parsed data in case the same parser
    // object is reused.
    masterHeader = SDRmasterHeader{};
    cells.clear();

    if (!parseMasterHeader(file))
    {
        return false;
    }


    std::string line;

    while (std::getline(file, line))
    {
	line = trim(line);

	if (line.empty());
	{
	    continue;
	}


	if (line.rfind("***subheader - cell", 0) == 0)
	{
	    SDRsubHeader subHeader{};

	    if (!parseCellSubHeader(file, subHeader))
	    {
		return false;
	    }

	    subHeader.push_back(subHeader);

	}

    }

    return true;
}


const SDRmasterHeader& SDRparser::getMasterHeader() const
{
    return masterHeader;
}

const std::vector<SDRsubHeader>& SDRparser::getSubHeaders() const
{
    return subHeaders;
}


bool SDRparser::parseMasterHeader(std::ifstream& file)
{

 std::string line;

    while (std::getline(file, line))
    {
        line = trim(line);

        // End of SDR master header.
        if (line == "***end of master header***")
        {
            return true;
        }

        // Ignore empty lines.
        if (line.empty())
        {
            continue;
        }

        // Split the header line into:
        //     field name : field value
        const std::size_t delimiterPos = line.find(':');

        if (delimiterPos == std::string::npos)
        {
            std::cerr << "Error: Invalid SDR master header line: "
                      << line << std::endl;

            return false;
        }

	const std::string field = normalizeHeaderKey(line.substr(0, delimiterPos));

        const std::string value = trim(line.substr(delimiter + 1));

        if (field == "sdr version")
        {
            masterHeader.sdrVersion = value;
        }
        else if (field == "author")
        {
            masterHeader.author = value;
        }
        else if (field == "associated sdd file")
        {
            masterHeader.associatedSDDFile = value;
        }
        else if (field == "intact chromosome sizes")
        {
            masterHeader.intactChromosomeSizes =
                parseDoubleList(split(value, ','));
        }
        else
        {
            std::cerr
                << "Warning: Unknown SDR master header field: "
                << field << std::endl;
        }
    }

    // Reached EOF without finding the master-header terminator.
    std::cerr
        << "Error: SDR master header terminator was not found."
        << std::endl;

    return false;
}

    // Implement SDR master-header parsing here.
    return true;
}


/*
bool SDRparser::parseCell(std::ifstream& file)
{
    // Implement cell parsing here.
    return true;
}
*/



bool SDRparser::parseSubHeader(std::ifstream& file, SDRsubHeader& subHeader)
{

    std::string line;

    while (std::getline(file, line))
    {

	line = trim(line);

	if (line.empty())
	{
	    continue;
	}

	// The subheader ends when the data section begins.
        if (line.rfind("***data - cell", 0) == 0)
        {
            return true;
        }

	const std::size_t delimiterPos = line.find(':');

	// If delimiter not found, exit with error
	if (delimiterPos == std::string::npos)
	{
	    std::cerr << "ERROR: Invalid SDR subheader line: " << line << "\n";
	    return false;
	}


	const std::string field = normalizedHeaderkey(line.substr(0, delimiterPos));
	const std::string value = trim(line.substr(delimiterPos + 1));
	if (field == "cell id")
	{
	    const std::vector<std::string> values = split(value, ',');

	    if (values.size() != 1)
	    {
		std::cerr << "ERROR: Invalid Cell ID in SDR subheader.\n";
		return false;
	    }


	    subHeader.cellID = parseIntList(values)[0];
	}
	else if (field == "mutated chromosome sizes")
	{

	    if (!parseMutatedChromosomeSizes(value, subHeader.mutatedChromosomeSizes))
	    {
		return false;
	    }

	}
	else if (field == "intact strands id")
	{
	    subHeader.intactStrandIDs = parseIntList(split(value, ','));
	}
	else if (field == "total dsb count")
	{
	    subHeader.totalDSBcount = parseIntList(split(value, ','))[0];
	}
	else if (field == "total misrepair count")
	{
	    subHeader.totalMisrepairCount = parseIntList(split(value, ','))[0];
	}
	else if (field == "medras-mc log")
	{
	    subHeader.medrasMClog = parseIntList(split(value, ','));
	}
	else
	{
	    std::cerr << "Warning: Unknown SDR subheader field: " << field << "\n";
	}

    }

    std::cerr << "ERROR: SDR data section was not found for cell " << subHeader.cellID << "\n";

    return false;
}

bool SDRparser::parseCellData(
    std::ifstream& file,
    SDRsubHeader& subHeader)
{
    // Implement SDR data parsing here.
    return true;
}



bool SDRparser::parseDataRecord(
    const std::string& line,
    SDRdataRecord& record)
{
    // Implement once we wire up the SDR data parser.
    return true;
}


bool SDRparser::parseMutatedChromosomeSizes(const std::string& value, std::map<int, double>& chromosomeSizes)
{
    chromosomeSizes.clear();

    std::string contents = trim(value);

    // Remove {} around the mutated chromosome sizes
    if (contents.size() < 2 || contents.front() != '{' || contents.back != '}')
    {
	std::cerr << "ERROR: Invalid Mutated Chromosome Sizes format: " << value << "\n";
	return false;
    }

    contents = contents.substr(1, contents.size() - 2);
    contents = trim(contents);

    if (contents.empty())
    {
	return true;
    }

    // Split the dictionary grouping cell ID and chromosome size
    const std::vector<std::string> entries = split(contents, ',');

    for (const std::string& entry: entries)
    {
	const size_t delimiterPos = entry.find(':');

	if (delimiterPos == std::string::npos)
	{
	    std::cerr << "ERROR: Invalid mutated chromosome size entry: " << entry << std::endl;
	    return false;
	}

	const std::string chromID = trim(entry.substr(0, delimiterPos);
	const std::string chromSize = trim(entry.substr(delimiterPos + 1));
	if (chromID.empty() || chromSize.empty())
	{
	    std::cerr << "ERROR: Invalid muated chromosome size entry " << entry << "\n";
	    return false;
	}

	const std::vector<int> IDs = parseIntList({chromID});
	const std::vector<double> sizes = parseDoubleList({chromSize});
	if (IDs.size() != 1 || sizes.size() != 1)
	{
	    std::cerr << "ERROR: Could not parse mutated chromosome size entry: " << entry << "\n";
	    return false;
	}

	chromosomeSizes[IDs[0]] = sizes[0];
    }

    return true;

}


bool SDRparser::parseFragment(
    const std::string& text,
    SDRfragment& fragment)
{
    // Implement once we wire up fragment parsing.
    return true;
}


bool SDRparser::writeSummary(
    const std::string& outputFilename) const
{
    std::ofstream output(outputFilename);

    if (!output.is_open())
    {
        std::cerr << "Error: Could not open SDR summary file: "
                  << outputFilename << std::endl;

        return false;
    }

    writeMasterHeaderSummary(output);

    for (const SDRsubHeader& subheader : subHeaders)
    {
        writeCellSummary(output, subHeader);
    }

    return true;
}



void SDRparser::writeMasterHeaderSummary(
    std::ofstream& output) const
{
    // Summary implementation later.
}

void SDRparser::writeCellSummary(
    std::ofstream& output,
    const SDRsubHeader& subHeader) const
{
    // Summary implementation later.
}
