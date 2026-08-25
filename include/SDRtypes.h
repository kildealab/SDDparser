#ifndef SDRTYPES_H
#define SDRTYPES_H

#include <string>
#include <vector>
#include <map>

// Represents one fragment used to construct a new SDR strand. - Field 3
struct SDRfragment
{
    int oldStrandID;
    double oldStartPosition;
    double oldEndPosition;
    bool hasCentromere;
};

// Represents one strand in the SDR data section.
struct SDRdataRecord
{
    int cellID;
    int newStrandID;
    std::vector<SDRfragment> fragments;

    // false = ring/circular
    // true  = linear
    bool linear;
};

// Stores the SDR master header.
struct SDRmasterHeader
{
    std::string sdrVersion;
    std::string author;
    std::string associatedSDDFile;

    // First value in the SDR file is the number of chromosomes,
    // followed by the chromosome sizes.
    std::vector<double> intactChromosomeSizes;
};

// Stores the information contained in one SDR cell subheader.
struct SDRsubHeader
{
    int cellID;

    // Key = mutated/new strand ID
    // Value = mutated chromosome/strand size in Mbp
    std::map<int, double> mutatedChromosomeSizes;

    // IDs of strands that remain intact.
    std::vector<int> intactStrandIDs;

    int totalDSBCount;
    int totalMisrepairCount;

    // Raw MEDRAS-MC log
    std::vector<int> medrasMClog;

    // All SDR data records belonging to this cell.
    std::vector<SDRdataRecord> dataRecords;

};

#endif
