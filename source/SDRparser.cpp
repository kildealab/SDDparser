#include "SDRparser.h"

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

    // Cell parsing will be implemented once we finalize
    // the SDR section markers and parsing behavior.

    return true;
}


const SDRmasterHeader& SDRparser::getMasterHeader() const
{
    return masterHeader;
}

const std::vector<SDRcell>& SDRparser::getCells() const
{
    return cells;
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



bool SDRparser::parseCellSubheader(
    std::ifstream& file,
    SDRcell& cell)
{
    // Implement cell subheader parsing here.
    return true;
}

bool SDRparser::parseCellData(
    std::ifstream& file,
    SDRcell& cell)
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

    for (const SDRcell& cell : cells)
    {
        writeCellSummary(output, cell);
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
    const SDRcell& cell) const
{
    // Summary implementation later.
}
