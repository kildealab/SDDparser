#ifndef SDR_UTILITIES_H
#define SDR_UTILITIES_H

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include "SDRtypes.h"

// -------------------------------------------------- //
// --------- SDR MUTATION DETECTION HELPERS --------- //
// -------------------------------------------------- //

// INVARIANT: SDRfragment::oldStrandID always refers to one of the
// original intact chromosomes (46 for humans), never to a
// newly-created strand ID produced by an earlier rearrangement.
// This is what allows every detect<Mutation>() function below to group a
// cell's fragments by oldStrandID and treat that as "all the
// rearrangement activity affecting a given original chromosome",
// rather than needing to trace strand lineage across events.

// Mutated strands are only those with new strand IDs >= numOriginalStrands.
// EX: for the 46 human chromosomes with strand IDs (1-46), mutated strands
// begin with new Strand ID > 46. This will be consistent across SDR files.
inline bool isRearrangementCandidate(int newStrandID, int numOriginalStrands)
{
    return newStrandID >= numOriginalStrands;
}


// 0.001 Mbp tolerance for long deletions, 10^-5 Mbp tolerance for balanced inversion/translocations.
inline bool approxEqual(double a, double b, double tolerance)
{

    return std::fabs(a - b) <= tolerance;
}



// --------------------------------------------------------------------------- //
// Detects LONG DELETION events within a single cell.
// --------------------------------------------------------------------------- //

// A long deletion event corresponds to two new data entries, one with a new strand ID and two flanking fragments with the same old strand ID
// and the other with a new strand ID with the central excised/deleted fragment with the same old strand ID as the other new strand.
inline std::vector<SDRdeletionEvent> detectDeletions(const SDRsubHeader& subHeader, int numOriginalStrands)
{

    double delTolerance = 0.001;						// Tolerance for the amount of discrepancy between neighbouring fragment start and end locations for long deletions is 1000 bases = 0.001 Mbp.

    std::vector<SDRdeletionEvent> deletions;					// Vector to hold deletion event objects

    std::map<int, std::vector<std::pair<int, SDRfragment>>> groupOldStrand;	// Group every fragment in this cell by the old strand ID it references, remembering which new-strand record it came from.

    for (const SDRdataRecord& record : subHeader.dataRecords)			// Loop through each cell subheader/data section to detect number of long deletions
    {
	if (!record.linear)							// Long deletions must have linear fragments
	{
	    continue;
	}

	if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))	// Make sure we are skipping the intact strands for mutation detection
    	{
            continue;
    	}

        for (const SDRfragment& fragment : record.fragments)			// Loop through all the fragments in SDR data entry field 3, store at each old strand ID the new strand ID and the corresponding fragment.
        {
            groupOldStrand[fragment.oldStrandID].push_back(
                {record.newStrandID, fragment});
        }
    }

    for (auto& [oldStrandID, entries] : groupOldStrand)				// Loop through the list of {oldStrandID, entry} pairs in groupOldStrand vector
    {
        std::map<int, std::vector<SDRfragment>> groupNewStrand;			// Group the new strand fragments with their corresponding strand IDs, starting from strand ID >= numOriginalStrands

        for (const auto& [newStrandID, fragment] : entries)			// Store the newStrandID (data field 2) and its corresponding fragments (data field 3)
        {
            groupNewStrand[newStrandID].push_back(fragment);
        }

        if (groupNewStrand.size() != 2) 					// Simple long deletions result in exactly two new data entries: one with two flanking fragments, and the other with the single excised fragment.
        {
            continue;
        }

	// Initialize strand IDs and their corresponding fragments
        int remainingStrandID = -1;
        int excisedStrandID = -1;
        std::vector<SDRfragment>* flankingFragments = nullptr;
        std::vector<SDRfragment>* excisedFragmentVec = nullptr;

        for (auto& [newStrandID, fragments] : groupNewStrand)			// Loop through each pair in groupNewStrand vector
        {
            if (fragments.size() == 2)						// If the entry has two flanking fragments, 
            {
                remainingStrandID = newStrandID;				// Track the newStrandID of the flanking fragments 
                flankingFragments = &fragments;					// Associate the the flanking fragments to the tracked newStrandID
            }
            else if (fragments.size() == 1)					// Fragment size of 1 corresponds to the excised/deleted DNA fragment
            {
                excisedStrandID = newStrandID;					// Track the excised fragment newStrandID
                excisedFragmentVec = &fragments;				// Track the corresponding fragment to this newStrandID
            }
        }

        if (flankingFragments == nullptr || excisedFragmentVec == nullptr)
        {
            continue;
        }

        // Sort the flanking fragments by start position so we can find the gap between them regardless of file order.
        std::sort(flankingFragments->begin(), flankingFragments->end(),[](const SDRfragment& a, const SDRfragment& b)
	    {
                return a.oldStartPosition < b.oldStartPosition;
            });

        const SDRfragment& lowerFragment = (*flankingFragments)[0];		// Flanking fragment prior to the first double-strand break leading to a deletion
        const SDRfragment& upperFragment = (*flankingFragments)[1];		// Flanking fragment proceeding the second double-strand break leading to a deletion
        const SDRfragment& excisedFragment = (*excisedFragmentVec)[0];		// The excised/deleted fragment

        const double gapStart = lowerFragment.oldEndPosition;
        const double gapEnd = upperFragment.oldStartPosition;

        if (gapEnd <= gapStart)							// Check if a gap was left from the deletion event (left flanking fragment end position must be less than right flanking fragment start position)
        {
            continue;
        }

	// Use a tolerance of ~1000 bases for a long deletion event to be tracked, excised fragment start and end locations need not necessarily line up perfectly with the flanking fragments strand ends. 
        if (!approxEqual(excisedFragment.oldStartPosition, gapStart, delTolerance) || !approxEqual(excisedFragment.oldEndPosition, gapEnd, delTolerance))
        {
            continue;
        }

	// Store the locations of the deletion and the strands involved for later Karyogram plotting
        SDRdeletionEvent event{};
        event.oldStrandID = oldStrandID;
        event.deletionStart = gapStart;
        event.deletionEnd = gapEnd;
        event.remainingStrandID = remainingStrandID;
        event.excisedStrandID = excisedStrandID;

        deletions.push_back(event);						// Store necessary deletion geometric information in the deletion vector.
    }

    return deletions;								// Return deletion information, for now summarize number of deletions with deletion.size() in the writeCellDataSummary function.
}













// Helper function to check for candidate balanced inversion fragments, where the start and end positions of a fragment are flipped.
inline bool isReversedFragment(const SDRfragment& fragment)
{
    return fragment.oldStartPosition > fragment.oldEndPosition;
}


// ---------------------------------------------------------------------------- //
// Detects BALANCED INVERSIONS within a single cell.
// ---------------------------------------------------------------------------- //

// A balanced inversion is recognized within a single new-strand record: fragments referencing the same old strand ID must form an
// intact region with no gaps (balanced) and exactly one of those fragments must be reversed.
inline std::vector<SDRinversionEvent> detectInversions(const SDRsubHeader& subHeader, int numOriginalStrands)
{

    double balInvTolerance = 0.00001;						// Tolerance for the amount of discrepancy between neighbouring fragment start and end locations is ~10 bases for balanced inversions = 10^-5 Mbp.

    std::vector<SDRinversionEvent> inversions;					// Vector to store Inversions events (strand IDs, strand start and end locations).

    for (const SDRdataRecord& record : subHeader.dataRecords)			// Loop through each line in the SDR data entries for a given cell.
    {
	if (!record.linear)							// Inversions must be composed of linear fragments
	{
	    continue;
	}

	if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))  // Make sure we are skipping the intact strands for mutation detection
    	{
            continue;
    	}

	std::map<int, std::vector<SDRfragment>> groupOldStrand;			// Vector to store old strand ID and corresponding fragments

        for (const SDRfragment& fragment : record.fragments)			// Loop through each fragment in a given SDR data entry/record
        {
            groupOldStrand[fragment.oldStrandID].push_back(fragment);		// Store the fragment corresponding to a given oldStrandID
        }

        for (auto& [oldStrandID, fragments] : groupOldStrand)			// Loop through each pair in the groupOldStrand vector.
        {
            if (fragments.size() < 2)						// Need at least two fragments referencing the same old strand for there to be any rearrangement to detect.
            {
                continue;
            }

            int reversedCount = 0;						// Use to make sure balanced inversion has exactly one reverse fragment, otherwise it is a different mutation.
            for (const SDRfragment& fragment : fragments)			// Check all fragments in a given SDR data entry/record.
            {
                if (isReversedFragment(fragment))				// Use helper function to check if the fragment start position > the end position.
                {
                    ++reversedCount;						// Increment number of inversions if detected.
                }
            }

            if (reversedCount != 1)						// A balanced inversion has exactly one reversed fragment.
            {
                continue;
            }

            struct NormalizedFragment						// Normalize fragment position to account for potential old strand start positions being greater than old strand end positions
            {
                double low;							// The lower fragment start position (closer to the start of the p arm)
                double high;							// The higher fragment end position (closer to the end of the q arm)
                bool reversed;							// Boolean to flag if the fragment is truly reverse/inverted.
            };

            std::vector<NormalizedFragment> normalized;				// Initialize normalized vector holding NormalizedFragment struct
            normalized.reserve(fragments.size());				// Reserve memory to store the size of the fragments

            for (const SDRfragment& fragment : fragments)			// Loop through each fragment in a data record/entry.
            {
                normalized.push_back({						// Store in the normalized position the locations of the fragment start and end positions, and whether it is reversed using the isReversedFragment helper function.
                    std::min(fragment.oldStartPosition, fragment.oldEndPosition),
                    std::max(fragment.oldStartPosition, fragment.oldEndPosition),
                    isReversedFragment(fragment)
                });
            }

	    // Sort all normalized fragments in each old Strand ID group pair according to their old strand start positions.
            std::sort(normalized.begin(), normalized.end(), [](const NormalizedFragment& a, const NormalizedFragment& b)
                {
                    return a.low < b.low;
                });


            bool contiguous = true;						// Minimal loss of bases (<= 10) can occur in balanced inversion, contiguous must be true
            for (std::size_t i = 0; i + 1 < normalized.size(); ++i)		// Loop through all normalized fragments per old Strand ID
            {
                if (!approxEqual(normalized[i].high, normalized[i + 1].low, balInvTolerance))
                {
                    contiguous = false;
                    break;
                }
            }

            if (!contiguous)
            {
                continue;
            }

	    // Find reversed fragment in normalized fragment
            const auto reversedIt = std::find_if(normalized.begin(), normalized.end(), [](const NormalizedFragment& f)
		{
		    return f.reversed;
		});

            if (reversedIt == normalized.end())						// If no reversed fragments found, skip to next data entry/record, do not record an inversion event
            {
                continue;
            }

	    // Record inversion event details for karyogram plotting later
            SDRinversionEvent event{};
            event.oldStrandID = oldStrandID;
            event.inversionStart = reversedIt->low;
            event.inversionEnd = reversedIt->high;
            event.newStrandID = record.newStrandID;

            inversions.push_back(event);						// Append inversion event details to the inversions vector
        }
    }

    return inversions;									// In writeCellDataSummary, print inversions.size() for the number of balanced inversions present in the SDR file.
}











// A candidate "half" of a balanced translocation: a record with exactly two fragments, referencing two distinct old strand IDs,
// both fragments in their original (non-reversed) orientation.
inline bool getTranslocationCandidateFragments(const SDRdataRecord& record, int numOriginalStrands, SDRfragment& fragmentA, SDRfragment& fragmentB)
{

    if (!record.linear)									// Fragments must be linear in balanced translocations
    {
	return false;
    }

    if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))		// Make sure we are skipping the intact strands for mutation detection
    {
	return false;
    }

    if (record.fragments.size() != 2)							// Each resulting new data record must contain exactly two fragments each, otherwise it is not a balanced translocation.
    {
        return false;
    }

    const SDRfragment& first = record.fragments[0];					// Track first balanced translocation fragment
    const SDRfragment& second = record.fragments[1];					// Track second balanced translocation fragment

    if (first.oldStrandID == second.oldStrandID)					// Balanced translocations occur on two different intact/old strand IDs, cannot have them on the same strand ID.
    {
        return false;
    }

    if (isReversedFragment(first) || isReversedFragment(second))			// Check if fragments are inverted, if so, not a balanced translocation mutation
    {
        return false;
    }

    fragmentA = first;
    fragmentB = second;

    return true;
}











// Checks whether two fragments referencing the same old strand ID meet exactly at a shared boundary within tolerance - i.e. they
// are the two pieces resulting from a single break in that strand. If so, populates breakpoint with the position of that break.
inline bool fragmentsShareBreakpoint(const SDRfragment& fragment1, const SDRfragment& fragment2, double& breakpoint)
{
    double balTransTolerance = 0.00001; 						// Balanced translocation mutations can have a loss of up to <= 10 bases = 10^-5 Mbp.

    if (fragment1.oldStrandID != fragment2.oldStrandID)					// The two fragments that are involved in balanced translocations cannot originate from the same original strand ID.
    {
        return false;
    }

    if (isReversedFragment(fragment1) || isReversedFragment(fragment2))			// Balanced translocation fragments cannot be reversed.
    {
        return false;
    }

    if (approxEqual(fragment1.oldEndPosition, fragment2.oldStartPosition, balTransTolerance)) // If fragment start and ends are within tolerance, record the location of the break for the first strands.
    {
        breakpoint = fragment1.oldEndPosition;
        return true;
    }

    if (approxEqual(fragment2.oldEndPosition, fragment1.oldStartPosition, balTransTolerance)) // If fragment start and ends are within tolerance, record the location of the break for the second strand.
    {
        breakpoint = fragment2.oldEndPosition;
        return true;
    }

    return false;
}













// ------------------------------------------------------------------------------------	//
// Function to detect BALANCED TRANSLOCATION entries in SDR file
// ------------------------------------------------------------------------------------ //

// Looks for a pair of new-strand records that, between them, reference exactly the same two old strand IDs (one fragment per
// old strand ID per record), where each old strand's two fragments meet cleanly at a single breakpoint within tolerance.
inline std::vector<SDRtranslocationEvent> detectTranslocations(const SDRsubHeader& subHeader, int numOriginalStrands)
{
    std::vector<SDRtranslocationEvent> translocations;					// Translocations vector to store translocation events

    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;			// Records vector to store the SDR data records/entries

    for (std::size_t i = 0; i < records.size(); ++i)					// Loop through the SDR data record/entries.
    {
        SDRfragment fragmentA1{};
        SDRfragment fragmentB1{};

        if (!getTranslocationCandidateFragments(records[i], numOriginalStrands, fragmentA1, fragmentB1)) // Check if the fragments in the data record are balanced translocation candidates using the helper function getTranslocationCandidateFragments
        {
            continue;
        }

        for (std::size_t j = i + 1; j < records.size(); ++j)				// Check the fragments in the data record directly after the current data record/entry to check for candidate translocation fragment.
        {
            SDRfragment fragmentA2{};
            SDRfragment fragmentB2{};

            if (!getTranslocationCandidateFragments(records[j], numOriginalStrands, fragmentA2, fragmentB2)) // Check if the fragments in the next following data record row are potential balanced translocation candidates.
            {
                continue;
            }

            // If both SDR data records i and j are translocation candidates, match up which fragment in record j shares an 
	    // old strand ID with which fragment in record i.
            SDRfragment matchA{};
            SDRfragment matchB{};

            if (fragmentA2.oldStrandID == fragmentA1.oldStrandID && fragmentB2.oldStrandID == fragmentB1.oldStrandID)
            {
                matchA = fragmentA2;
                matchB = fragmentB2;
            }
            else if (fragmentA2.oldStrandID == fragmentB1.oldStrandID && fragmentB2.oldStrandID == fragmentA1.oldStrandID)
            {
                matchA = fragmentB2;
                matchB = fragmentA2;
            }
            else									// If old strand IDs do not match, not a balanced translocation, move on to next entry.
            {
                continue;
            }

	    // Track balanced translocation breakpoint locations in each data record
            double breakpointA = 0.0;
            double breakpointB = 0.0;

	    // If the breakpoints do not match up within tolerance (10 bases), they are not considered balanced translocations.
            if (!fragmentsShareBreakpoint(fragmentA1, matchA, breakpointA) || !fragmentsShareBreakpoint(fragmentB1, matchB, breakpointB))
            {
                continue;
            }

	    // All checks have passes, record the translocation event for future karyogram plotting
            SDRtranslocationEvent event{};
            event.oldStrandA = fragmentA1.oldStrandID;
            event.oldStrandB = fragmentB1.oldStrandID;
            event.breakpointA = breakpointA;
            event.breakpointB = breakpointB;
            event.newStrandID1 = records[i].newStrandID;
            event.newStrandID2 = records[j].newStrandID;

            translocations.push_back(event);						// Store the balanced translocation event in the translocations vector
        }
    }

    return translocations;								// In writeCellDataSummary, return the number of translocation events detected in the SDR file
}









// ----------------------------------------------------------------- //
// Detects EXTRACHROMOSOMAL DNA (ecDNA) events within a single cell. //
// ----------------------------------------------------------------- //

// Identical fragment signature to the long deletion events with the added
// indicator of data field 4 isLinear = 0 (for circular fragments). one new strand
// retains the two flanking fragments, another new strand is the
// excised middle segment - except the excised strand's record must
// be CIRCULAR (linear = 0) rather than linear. The remaining/flanking
// strand is still required to be linear.

inline std::vector<SDRecDNAevent> detectECDNA(const SDRsubHeader& subHeader, int numOriginalStrands)
{
    std::vector<SDRecDNAevent> ecDNAevents;

    // Group every fragment in this cell by the old strand ID it
    // references, remembering which new-strand record it came from.
    // Unlike detectDeletions, we can't filter out circular records
    // up front - the excised piece here is EXPECTED to be circular.
    std::map<int, std::vector<std::pair<const SDRdataRecord*, SDRfragment>>> groupOldStrand;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
	if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))
	{
	    continue;
	}

	for (const SDRfragment& fragment : record.fragments)
	{
	    groupOldStrand[fragment.oldStrandID].push_back({&record, fragment});
	}

    }

    for (auto& [oldStrandID, entries] : groupOldStrand)
    {
        std::map<int, std::vector<SDRfragment>> groupNewStrand;
        std::map<int, const SDRdataRecord*> recordByNewStrand;

        for (const auto& [record, fragment] : entries)
        {
            groupNewStrand[record->newStrandID].push_back(fragment);
            recordByNewStrand[record->newStrandID] = record;
        }

        // An ecDNA event touches exactly two new-strand records: one
        // linear record with the two flanking fragments, one
        // CIRCULAR record with the single excised fragment.
        if (groupNewStrand.size() != 2)
        {
            continue;
        }

	int remainingStrandID = -1;
        int excisedStrandID = -1;
        std::vector<SDRfragment>* flankingFragments = nullptr;
        std::vector<SDRfragment>* excisedFragmentVec = nullptr;

        for (auto& [newStrandID, fragments] : groupNewStrand)
        {
            const SDRdataRecord* record = recordByNewStrand[newStrandID];

            if (fragments.size() == 2 && record->linear)
            {
                remainingStrandID = newStrandID;
                flankingFragments = &fragments;
            }
            else if (fragments.size() == 1 && !record->linear)
            {
                excisedStrandID = newStrandID;
                excisedFragmentVec = &fragments;
            }
        }


	if (flankingFragments == nullptr || excisedFragmentVec == nullptr)
        {
            continue;
        }

        std::sort(flankingFragments->begin(), flankingFragments->end(), [](const SDRfragment& a, const SDRfragment& b)
        {
            return a.oldStartPosition < b.oldStartPosition;
        });

        const SDRfragment& lowerFragment = (*flankingFragments)[0];
        const SDRfragment& upperFragment = (*flankingFragments)[1];
        const SDRfragment& excisedFragment = (*excisedFragmentVec)[0];

	const double ecDNAtolerance = 0.001;					// tolerance in Mbp difference between contiguous segments <= 100 bp

        const double gapStart = lowerFragment.oldEndPosition;
        const double gapEnd = upperFragment.oldStartPosition;

	if (gapEnd <= gapStart)
        {
            continue;
        }

        if (!approxEqual(excisedFragment.oldStartPosition, gapStart, ecDNAtolerance) ||
            !approxEqual(excisedFragment.oldEndPosition, gapEnd, ecDNAtolerance))
        {
            continue;
        }

	if (excisedFragment.hasCentromere)					// ecDNA fragments are acentric (excised fragment should not contain a centromere)
	{
	    continue;
	}


        SDRecDNAevent event{};
        event.oldStrandID = oldStrandID;
        event.ecDNAstart = gapStart;
        event.ecDNAend = gapEnd;
        event.remainingStrandID = remainingStrandID;
        event.excisedStrandID = excisedStrandID;

        ecDNAevents.push_back(event);
    }

    return ecDNAevents;

}












// ---------------------------------------------------------------- //
// Function to detect DELETION-INVERSION mutations in a single cell //
// ---------------------------------------------------------------- //

// One new-strand record contains three fragments from the same old
// strand: two in normal orientation and one reversed. A second
// new-strand record contains a single fragment from the same old
// strand - the deleted segment. Unlike a plain balanced inversion,
// the three fragments are NOT fully contiguous: there is exactly one
// gap among them (where the deletion occurred), and that gap must
// exactly match the second record's excised fragment (within tolerance).
inline std::vector<SDRdeletionInversionEvent> detectDeletionInversions(
    const SDRsubHeader& subHeader,
    int numOriginalStrands)
{
    std::vector<SDRdeletionInversionEvent> delInvs;

    std::map<int, std::vector<std::pair<int, SDRfragment>>> groupOldStrand;

    const double delInvTolerance = 0.001;					// ~ 0.001 Mbp = 1000 bp tolerance for position mismatch

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (!record.linear)
        {
            continue;
        }

        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))
        {
            continue;
        }

        for (const SDRfragment& fragment : record.fragments)
        {
            groupOldStrand[fragment.oldStrandID].push_back({record.newStrandID, fragment});
        }
    }

    for (auto& [oldStrandID, entries] : groupOldStrand)
    {
        std::map<int, std::vector<SDRfragment>> groupNewStrand;

        for (const auto& [newStrandID, fragment] : entries)
        {
            groupNewStrand[newStrandID].push_back(fragment);
        }

        // Must touch exactly two new-strand records: one with three
        // fragments (2 normal + 1 inverted), one with the single
        // excised fragment.
        if (groupNewStrand.size() != 2)
        {
            continue;
        }

        int remainingStrandID = -1;
        int excisedStrandID = -1;
        std::vector<SDRfragment>* threeFragments = nullptr;				// Check if one new strand has 3 fragments
        std::vector<SDRfragment>* excisedFragmentVec = nullptr;				// Check if one new strand has 1 fragment (deleted portion)

        for (auto& [newStrandID, fragments] : groupNewStrand)
        {
            if (fragments.size() == 3)							// If this entry has 3 fragments, store it in the remainingStrandID vector (has two original strands and one inversion)
            {
                remainingStrandID = newStrandID;
                threeFragments = &fragments;
            }
            else if (fragments.size() == 1)						// If this entry has 1 fragment, store it in the excisedFragmentID vector. 
            {
                excisedStrandID = newStrandID;
                excisedFragmentVec = &fragments;
            }
        }

        if (threeFragments == nullptr || excisedFragmentVec == nullptr)			// If the 3 fragment - 1 fragment structure not observed, not a deletion-inversion
        {
            continue;
        }

        const SDRfragment& excisedFragment = (*excisedFragmentVec)[0];

        if (isReversedFragment(excisedFragment))
        {
            continue; 									// The excised piece itself should not be reversed.
        }


        int reversedCount = 0;

        for (const SDRfragment& fragment : *threeFragments)				// Check how many inversions in the remainingStrandID strand
        {
            if (isReversedFragment(fragment))
            {
                ++reversedCount;
            }
        }

        if (reversedCount != 1)								// If more than one inversion, not considered a deletion-inversion event.
        {
            continue;
        }

        struct NormalizedFragment							// Local structure to order the fragments
        {
            double low;					// Stores the lower of the two positions of the fragment (not necessarily the start position in an inverted fragment)
            double high;				// Stores the higher of the two positions of the fragment (not necessarily the end position in an inversion)
            bool reversed;
        };

        std::vector<NormalizedFragment> normalizedFragment;
        normalizedFragment.reserve(3);								// Reserve at least 3 elements for normalized fragment

	// Reorder all fragments in terms of start and end locations (inversions become uninverted to sort properly)
        for (const SDRfragment& fragment : *threeFragments)
        {
            normalizedFragment.push_back({
                std::min(fragment.oldStartPosition, fragment.oldEndPosition),
                std::max(fragment.oldStartPosition, fragment.oldEndPosition),
                isReversedFragment(fragment)
            });
        }

	// Sort all fragments in order of their start positions in ascending order.
        std::sort(normalizedFragment.begin(), normalizedFragment.end(), [](const NormalizedFragment& a, const NormalizedFragment& b)
        {
            return a.low < b.low;
        });

        // Check the two boundaries between the three sorted fragments -> exactly one should have a gap (the deletion), the other
        // should be contiguous.
        double gapStart = 0.0;						// In the record with the three fragments, record where the deletion started and ended
        double gapEnd = 0.0;
        int gapCount = 0;						// Count number of gaps
        bool overlapFound = false;					// Check if gap overlaps with a fragment in another new strand record (where the excised fragment is)
        for (std::size_t i = 0; i + 1 < normalizedFragment.size(); ++i)
        {
            if (approxEqual(normalizedFragment[i].high, normalizedFragment[i + 1].low, delInvTolerance))
            {
                continue; // Contiguous - no gap here.
            }

            if (normalizedFragment[i + 1].low <= normalizedFragment[i].high)
            {
                overlapFound = true; // Overlapping, not a gap - invalid shape.
                break;
            }

            gapStart = normalizedFragment[i].high;
            gapEnd = normalizedFragment[i + 1].low;
            ++gapCount;
        }

        if (overlapFound || gapCount != 1)
        {
            continue;
        }

        if (!approxEqual(excisedFragment.oldStartPosition, gapStart, delInvTolerance) ||
            !approxEqual(excisedFragment.oldEndPosition, gapEnd, delInvTolerance))
        {
            continue;
        }

        const auto reversedIt = std::find_if(normalizedFragment.begin(), normalizedFragment.end(), [](const NormalizedFragment& f)
        {
            return f.reversed;
        });

        if (reversedIt == normalizedFragment.end())
        {
            continue;
        }

        SDRdeletionInversionEvent event{};
        event.oldStrandID = oldStrandID;
        event.inversionStart = reversedIt->low;
        event.inversionEnd = reversedIt->high;
        event.deletionStart = gapStart;
        event.deletionEnd = gapEnd;
        event.remainingStrandID = remainingStrandID;
        event.excisedStrandID = excisedStrandID;

        delInvs.push_back(event);
    }

    return delInvs;
}











// Checks whether two fragments referencing the same old strand ID
// have a genuine GAP between them (rather than meeting cleanly at a
// breakpoint) - i.e. material between them is missing. If so,
// populates gapStart/gapEnd with the gap's bounds.
inline bool fragmentsHaveGap(const SDRfragment& fragment1, const SDRfragment& fragment2, double& gapStart, double& gapEnd)
{
    const double delTolerance = 0.001; 				// Tolerance for base pair positioning mismatch for deletions is 0.001 Mbp = 1000 bp or less

    if (fragment1.oldStrandID != fragment2.oldStrandID)
    {
        return false;
    }

    if (isReversedFragment(fragment1) || isReversedFragment(fragment2))
    {
        return false;
    }

    if (fragment1.oldEndPosition < fragment2.oldStartPosition &&
        !approxEqual(fragment1.oldEndPosition, fragment2.oldStartPosition, delTolerance))
    {
        gapStart = fragment1.oldEndPosition;
        gapEnd = fragment2.oldStartPosition;
        return true;
    }

    if (fragment2.oldEndPosition < fragment1.oldStartPosition &&
        !approxEqual(fragment2.oldEndPosition, fragment1.oldStartPosition, delTolerance))
    {
        gapStart = fragment2.oldEndPosition;
        gapEnd = fragment1.oldStartPosition;
        return true;
    }

    return false;
}











// --------------------------------------------------------------- //
// Function to detect BALANCED TRANSLOCATIONS in a single cell     //
// --------------------------------------------------------------- //

// Looks for 2 fragments each, from 2 different old strands, linear, 
// non-reversed. If BOTH strands meet
// cleanly at a breakpoint, that's a plain balanced translocation -
// not this. If exactly ONE strand meets cleanly and the other has a
// gap, a third record (single fragment, from the gapped strand,
// exactly filling the gap) confirms a deletion-translocation.

inline std::vector<SDRdeletionTranslocationEvent> detectDeletionTranslocations(
    const SDRsubHeader& subHeader,
    int numOriginalStrands)
{
    std::vector<SDRdeletionTranslocationEvent> delTras;

    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;

    const double delTraTolerance = 0.001;					// Tolerance for the mismathc in contiguous base positions is 0.001 Mbp = 1000 bp or less

    for (std::size_t i = 0; i < records.size(); ++i)
    {
        SDRfragment fragmentA1{};
        SDRfragment fragmentB1{};

        if (!getTranslocationCandidateFragments(records[i], numOriginalStrands, fragmentA1, fragmentB1))
        {
            continue;
        }

        for (std::size_t j = i + 1; j < records.size(); ++j)
        {
            SDRfragment fragmentA2{};
            SDRfragment fragmentB2{};

            if (!getTranslocationCandidateFragments(records[j], numOriginalStrands, fragmentA2, fragmentB2))
            {
                continue;
            }

            SDRfragment matchA{};
            SDRfragment matchB{};

            if (fragmentA2.oldStrandID == fragmentA1.oldStrandID && fragmentB2.oldStrandID == fragmentB1.oldStrandID)
            {
                matchA = fragmentA2;
                matchB = fragmentB2;
            }
            else if (fragmentA2.oldStrandID == fragmentB1.oldStrandID && fragmentB2.oldStrandID == fragmentA1.oldStrandID)
            {
                matchA = fragmentB2;
                matchB = fragmentA2;
            }
            else
            {
                continue;
            }

            double breakpointA = 0.0;
            const bool cleanA = fragmentsShareBreakpoint(fragmentA1, matchA, breakpointA);

            double breakpointB = 0.0;
            const bool cleanB = fragmentsShareBreakpoint(fragmentB1, matchB, breakpointB);

            if (cleanA && cleanB)
            {
                continue; // Plain balanced translocation - handled elsewhere, not here.
            }

            double gapStart = 0.0;
            double gapEnd = 0.0;
            int deletedOldStrandID = -1;
            double cleanBreakPos = 0.0;

            if (cleanA && !cleanB)
            {
                if (!fragmentsHaveGap(fragmentB1, matchB, gapStart, gapEnd))
                {
                    continue;
                }

                deletedOldStrandID = fragmentB1.oldStrandID;
                cleanBreakPos = breakpointA;
            }
            else if (cleanB && !cleanA)
            {
                if (!fragmentsHaveGap(fragmentA1, matchA, gapStart, gapEnd))
                {
                    continue;
                }

                deletedOldStrandID = fragmentA1.oldStrandID;
                cleanBreakPos = breakpointB;
            }
            else
            {
                continue; // Neither strand cleanly matches - not this shape.
            }

            for (const SDRdataRecord& thirdRecord : records)
            {
                if (thirdRecord.newStrandID == records[i].newStrandID || thirdRecord.newStrandID == records[j].newStrandID)
                {
                    continue;
                }

                if (!thirdRecord.linear)
                {
                    continue;
                }

                if (!isRearrangementCandidate(thirdRecord.newStrandID, numOriginalStrands))
                {
                    continue;
                }

                if (thirdRecord.fragments.size() != 1)
                {
                    continue;
                }

                const SDRfragment& excisedFragment = thirdRecord.fragments[0];

                if (excisedFragment.oldStrandID != deletedOldStrandID)
                {
                    continue;
                }

                if (isReversedFragment(excisedFragment))
                {
                    continue;
                }

                if (!approxEqual(excisedFragment.oldStartPosition, gapStart, delTraTolerance) ||
                    !approxEqual(excisedFragment.oldEndPosition, gapEnd, delTraTolerance))
                {
                    continue;
                }

                SDRdeletionTranslocationEvent event{};
                event.oldStrandA = fragmentA1.oldStrandID;
                event.oldStrandB = fragmentB1.oldStrandID;
                event.deletedOldStrandID = deletedOldStrandID;
                event.cleanBreakPos = cleanBreakPos;
                event.deletionStart = gapStart;
                event.deletionEnd = gapEnd;
                event.newStrandID1 = records[i].newStrandID;
                event.newStrandID2 = records[j].newStrandID;
                event.excisedStrandID = thirdRecord.newStrandID;

                delTras.push_back(event);
                break;
            }
        }
    }

    return delTras;
}











#endif
