#ifndef SDR_TYPES_H
#define SDR_TYPES_H

#include <string>
#include <vector>
#include <map>

// Represents one fragment used to construct a new SDR strand. - Field 3
struct SDRfragment
{
// Field 3: oldStrandID/oldStrandStart/oldStrandEnd/hasCentromere
    int oldStrandID;
    double oldStartPosition;
    double oldEndPosition;
    bool hasCentromere;
};

// Represents one strand in the SDR data section.
struct SDRdataRecord
{
    int cellID;						// SDR data field 1
    int newStrandID;					// SDR data field 2
    std::vector<SDRfragment> fragments;			// SDR data field 3, can contain multiple subfields separated by ',' and each subfield contains 4 entries separated by '/'.

    // false = ring/circular
    // true  = linear
    bool linear;
};

// Stores the SDR master header entries
struct SDRmasterHeader
{
    std::string sdrVersion;
    std::string author;
    std::string associatedSDDFile;
    std::vector<double> intactChromosomeSizes;		// First value in the SDR file is the number of chromosomes, followed by the chromosome sizes.
};

// Stores the information contained in one SDR cell subheader.
struct SDRsubHeader
{
    int cellID;

    std::map<int, double> mutatedChromosomeSizes;	// Key = mutated/new strand ID, Value = mutated chromosome/strand size in Mbp

    std::vector<int> intactStrandIDs;			// IDs of strands that remain intact.

    int totalDSBcount;					// A blank value in SDR subheader means not measured, not an error
    int totalMisrepairCount;

    std::vector<int> medrasMClog;

    std::vector<SDRdataRecord> dataRecords;		// All SDR data records belonging to this cell.

};

// Use for checking if total DSB count, total Misrepair count, or MEDRAS-MC log are blank, if so pass warning, not error.
constexpr int SDR_FIELD_NOT_MEASURED = -1;


// ------------------------------------------------ //
// Mutation types for summarizing in text file     
// ------------------------------------------------ //


// Represents a detected deletion event within a single cell
struct SDRdeletionEvent
{
    int oldStrandID;
    double deletionStart;
    double deletionEnd;
    int remainingStrandID;   // new strand retaining the two flanking fragments
    int excisedStrandID;     // new strand consisting of just the deleted segment
};


// Represents a detected balanced inversion: a segment of
// the original strand was flipped in orientation in place, with no
// net gain or loss of material.
struct SDRinversionEvent
{
    int oldStrandID;
    double inversionStart;
    double inversionEnd;
    int newStrandID;
};


// Represents a detected balanced translocation: two
// original strands were each split once, and the resulting pieces
// were swapped between two new strands, with no material gained or
// lost.
struct SDRtranslocationEvent
{
    int oldStrandA;
    int oldStrandB;
    double breakpointA;   // position on oldStrandA where the break occurred
    double breakpointB;   // position on oldStrandB where the break occurred
    int newStrandID1;
    int newStrandID2;
};




#endif
