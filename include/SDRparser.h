#ifndef SDRPARSER_H
#define SDRPARSER_H

#include "SDRtypes.h"

#include <fstream>
#include <string>
#include <vector>


class SDRparser
{
public:

    SDRparser() = default;

    // Parse an SDR file.
    bool parseFile(const std::string& filename);

    // Access parsed SDR information.
    const SDRmasterHeader& getMasterHeader() const;
    const std::vector<SDRsubHeader>& getSubHeaders() const;

    // Generate a summary file.
    bool writeSummary(const std::string& outputFilename) const;

private:

    SDRmasterHeader masterHeader;
    std::vector<SDRsubHeader> subHeaders;

    // Master header parsing.
    bool parseMasterHeader(std::ifstream& file);
/*
    // Cell parsing.
    bool parseCell(std::ifstream& file);
*/
    // Cell subheader parsing.
    bool parseCellSubheader(
        std::ifstream& file,
        SDRsubHeader& subHeader
    );

    // Cell data parsing.
    bool parseCellData(
        std::ifstream& file,
        SDRsubHeader& subHeader
    );


    bool parseMutatedChromosomeSizes(const std::string& value, std::map<int, double>& chromosomeSizes);


    // Parse one SDR data record.
    bool parseDataRecord(
        const std::string& line,
        SDRdataRecord& record
    );

    // Parse one fragment from field 3.
    bool parseFragment(
        const std::string& text,
        SDRfragment& fragment
    );

    // Summary helpers.
    void writeMasterHeaderSummary(std::ofstream& output) const;
    void writeSubHeaderSummary(
        std::ofstream& output,
        const SDRsubHeader& subHeader
    ) const;
};

#endif
