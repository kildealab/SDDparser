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
    int oldStrandID;		// Which chromosome/strand ID the fragments originate from
    double deletionStart;	// Location (in Mbp from start of p arm) for the beginning of the deletion
    double deletionEnd;		// Location (in Mbp from start of p arm) for the beginning of the deletion
    int remainingStrandID;   	// new strand retaining the two flanking fragments
    int excisedStrandID;     	// new strand consisting of just the deleted segment
};


// Represents a detected balanced inversion: a segment of
// the original strand was flipped in orientation in place, with no
// net gain or loss of material.
struct SDRinversionEvent
{
    int oldStrandID;		// Which chromosome/strand ID the fragments originate from
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
    int oldStrandA;		// Which chromosome/strand ID the fragments originate from for strand A
    int oldStrandB;		// Which chromosome/strand ID the fragments originate from for strand B
    double breakpointA;   	// position on oldStrandAid where the break occurred
    double breakpointB;   	// position on oldStrandBid where the break occurred
    int newStrandID1;		// ID of the new strand containing fragment rearrangements for strand 1
    int newStrandID2;		// ID of the new strand containing fragment rearrangenemts for strand 2
};


// Represents a detected extrachromosomal DNA fragment -
// It has the same fragment signature as long deletions, but
// the isLinear data field 4 is 0 (for circular fragments) rather than
// 1 for linear fragments.
struct SDRecDNAevent
{
    int oldStrandID;		// Which chromosome/strand ID the fragments originate from
    double ecDNAstart;		// position on oldStrandID where the break started (in Mbp from the start of the p arm)
    double ecDNAend;		// position on oldStrandID where the break ended (in Mbp from the start of the p arm)
    int remainingStrandID;	// New strand ID of the strand that lost the excised fragment
    int excisedStrandID;	// New strand ID of the excised strand
};




// Represents a combined detected deletion and inversion event -
// One new strand contains three fragments with an inversion and a
// missing section. Another new strand contains one fragment that is
// the missing section from the previous new strand.
struct SDRdeletionInversionEvent
{
    int oldStrandID;		// Which chromosome/strand ID the fragments originate from
    double deletionStart;	// Location of the start of the deletion in Mbp
    double deletionEnd;		// Location of the end of the deletion in Mbp
    double inversionStart;	// Location of the start of the inversion in Mbp
    double inversionEnd;	// Location of the end of the inversion in Mbp
    int remainingStrandID;	// New strand containing the three fragments, 2 original and one inverted.
    int excisedStrandID;	// New strand containing the deleted section, one fragment.
};





// Represents a combined detected deletion and translocation event - 
// Two strands contain two translocated fragments, and another strand
// contains a deleted segment of one of the original strand IDs, for 
// a total of three new data records per delTra mutation
struct SDRdeletionTranslocationEvent
{
    int oldStrandA;		// Old strand ID for strand A involved in the translocation
    int oldStrandB;		// Old strand ID for strand B involved in the translocation
    int deletedOldStrandID;   	// whichever of oldStrandA/oldStrandB has the gap, the
    double cleanBreakPos;     	// breakpoint is on the OTHER strand, which split cleanly
    double deletionStart;	// Location of where the deletion started in Mbp
    double deletionEnd;		// Location of where the deletion ended in Mbp
    int newStrandID1;         	// First translocation record ID
    int newStrandID2;         	// Second translocation record ID
    int excisedStrandID;      	// The excised strand's ID

};

#endif
