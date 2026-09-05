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

// Modified detectDeletions function to detect multiple deletions within a single strand. A strand with N deletions will contain N + 1 fragments in SDR data field 3, and will have
// N + 1 data entries, 1 entry being the original strand with all the gaps, and N entries representing each deletion causing the gaps.
inline std::vector<SDRdeletionEvent> detectDeletions(const SDRsubHeader& subHeader, int numOriginalStrands)
{

    double delTolerance = 0.001;						// Tolerance for the amount of discrepancy between neighbouring fragment start and end locations for long deletions is 1000 bases = 0.001 Mbp.

    std::vector<SDRdeletionEvent> deletions;					// Vector to hold deletion event objects

    std::map<int, std::vector<std::pair<int, SDRfragment>>> groupOldStrand;	// Group every fragment in this cell by the old strand ID it references, remembering which new-strand record it came from.


    // groupOldStrand looks like this:

    // groupOldStrand[oldStrandID] = { { &recordnewStrandID, Fragment{oldStrandID, start, end, hasCentromere} }, ... (times N) }

    // where N is the number of rearranged fragments that reference a given oldStrandID (i.e. for a newStrandId > 46, if three fragments reference a given
    // oldStrandId, then groupOldStrand[oldStrandID] will have three entries resembling the above format.


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
            groupOldStrand[fragment.oldStrandID].push_back({record.newStrandID, fragment});
        }
    }


    for (auto& [oldStrandID, entries] : groupOldStrand)				// Loop through the list of {oldStrandID, entry} pairs in groupOldStrand vector
    {
        std::map<int, std::vector<SDRfragment>> groupNewStrand;			// Group the new strand fragments with their corresponding strand IDs, starting from strand ID >= numOriginalStrands

        for (const auto& [newStrandID, fragment] : entries)			// Store the newStrandID (data field 2) and its corresponding fragments (data field 3)
        {
            groupNewStrand[newStrandID].push_back(fragment);
        }

        if (groupNewStrand.size() < 2)  					// N deletions lead to N + 1 data entries, so if there is 1 deletion, there should be at least two data entries.
        {
            continue;
        }


	// groupNewStrand looks like this:

	// groupNewStrand[newStrandID] = { Fragment{start, end, hasCentromere}, ... (N times) }

	// N is the number of fragments contained in the data record of new strand ID.


	// Initialize strand IDs and their corresponding fragments
        int remainingStrandID = -1;
        std::vector<SDRfragment>* flankingFragments = nullptr;
	std::vector<std::pair<int, SDRfragment*>> excisedCandidates;

	bool multipleRemainingRecords = false;

        for (auto& [newStrandID, fragments] : groupNewStrand)			// Loop through each pair in groupNewStrand vector
        {
            if (fragments.size() >= 2)						// If the entry has two flanking fragments, 
            {
		if (flankingFragments != nullptr)
		{
		    multipleRemainingRecords = true;
		    break;
		}

		remainingStrandID = newStrandID;                                // Track the newStrandID of the flanking fragments
                flankingFragments = &fragments;
            }
            else if (fragments.size() == 1)					// Fragment size of 1 corresponds to the excised/deleted DNA fragment
            {
		excisedCandidates.push_back({newStrandID, &fragments[0]});

            }
        }

	if (multipleRemainingRecords || flankingFragments == nullptr || excisedCandidates.empty())
        {
            continue;
        }


	// N deletions produce N+1 surviving fragments and N excised pieces.
        if (flankingFragments->size() != excisedCandidates.size() + 1)
        {
            continue;
        }


        // Sort the flanking fragments by start position so we can find the gap between them regardless of file order.
        std::sort(flankingFragments->begin(), flankingFragments->end(),[](const SDRfragment& a, const SDRfragment& b)
	    {
                return a.oldStartPosition < b.oldStartPosition;
            });


	struct Gap
	{
	    double startPos;
	    double endPos;
	};

	std::vector<Gap> gaps;


	for (std::size_t i = 0; i + 1 < flankingFragments->size(); i++)
        {
            const double gapStart = (*flankingFragments)[i].oldEndPosition;
            const double gapEnd = (*flankingFragments)[i + 1].oldStartPosition;

            if (gapEnd <= gapStart)                                             // Check if a gap was left from the deletion event
            {
                gaps.clear();
                break;                                                          // Overlapping/non-gap boundary - not a valid deletion shape for this strand.
            }

            gaps.push_back({gapStart, gapEnd});
        }

	if (gaps.size() != excisedCandidates.size())                            // Gap count must match the number of excised records.
        {
            continue;
        }


	// Match each gap to the excised record whose fragment fills it
        // exactly (within tolerance), and build one event per match.
        std::vector<bool> excisedClaimed(excisedCandidates.size(), false);
        std::vector<SDRdeletionEvent> deletionEvents;
        bool allGapsMatched = true;

	for (const Gap& gap : gaps)
        {
            bool matched = false;

            for (std::size_t i = 0; i < excisedCandidates.size(); ++i)
            {
                if (excisedClaimed[i])
                {
                    continue;
                }

                const SDRfragment& excisedFragment = *excisedCandidates[i].second;

                // Use a tolerance of ~1000 bases for a long deletion event to be tracked.
                if (approxEqual(excisedFragment.oldStartPosition, gap.startPos, delTolerance) &&
                    approxEqual(excisedFragment.oldEndPosition, gap.endPos, delTolerance))
                {
                    excisedClaimed[i] = true;
                    matched = true;

                    // Store the locations of the deletion and the strands involved for later Karyogram plotting
                    SDRdeletionEvent event{};
		    event.oldStrandID = oldStrandID;
                    event.deletionStart = gap.startPos;
                    event.deletionEnd = gap.endPos;
                    event.remainingStrandID = remainingStrandID;
                    event.excisedStrandID = excisedCandidates[i].first;

                    deletionEvents.push_back(event);
                    break;
                }
            }

            if (!matched)
            {
                allGapsMatched = false;
                break;
            }
        }

	if (!allGapsMatched)                                                    // Some gap didn't match any excised record - not a valid deletion shape.
        {
            continue;
        }

        deletions.insert(deletions.end(), deletionEvents.begin(), deletionEvents.end());
    }

    return deletions;  

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
// indicator of data field 4 isLinear = 0 (for circular fragments). 
// Generalized to N independent ecDNA mutations, and for ecDNA made up of
// multiple deleted fragments from the same original strand. Can have either
// One remaining strand with N deletions and N+1 fragments, followed by
// N records of ecDNA depicting the N ecDNA mutations, or followed by one or 
// multiple records with multiple fragments forming ecDNA.

inline std::vector<SDRecDNAevent> detectECDNA(const SDRsubHeader& subHeader, int numOriginalStrands)
{
    std::vector<SDRecDNAevent> ecDNAevents;					// Vector to store ecDNA event information (strand IDs, fragment lengths)

    // Group every fragment in this cell by the old strand ID it
    // references, remembering which new-strand record it came from.
    // Unlike detectDeletions, the excised piece here is EXPECTED to be circular.
    std::map<int, std::vector<std::pair<const SDRdataRecord*, SDRfragment>>> groupOldStrand;

    for (const SDRdataRecord& record : subHeader.dataRecords)			// Loop through each data record per cell subheader
    {
	if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))	// Check if a data record consists of a rearranged fragment, if not, skip
	{
	    continue;
	}

	for (const SDRfragment& fragment : record.fragments)			// Check the fragments the new strand is composed of
	{
	    groupOldStrand[fragment.oldStrandID].push_back({&record, fragment});// Store the fragments corresponding to each original (old) strand ID
	}
    }


    // This regroups the same three fragments, but now keyed by newStrandID instead of oldStrandID,
    // collapsing "which record" down to just an ID, and separately keeping a lookup back to the actual record object:
    for (auto& [oldStrandID, entries] : groupOldStrand)				// Loop through all the entries/fragments associated with a given old strand ID
    {
        std::map<int, std::vector<SDRfragment>> groupNewStrand;			// Vector to store the new strand ID and the fragments the new strand is composed of -> use to tell how many fragments each new strand ID is made of
        std::map<int, const SDRdataRecord*> recordByNewStrand;			// Need to get back from newStrandID the actual SDR data record to check if strand is linear or circular

        for (const auto& [record, fragment] : entries)
        {
            groupNewStrand[record->newStrandID].push_back(fragment);
            recordByNewStrand[record->newStrandID] = record;
        }


	// recordByNewStrand looks like this:

	// recordByNewStrand[newStrandID] = &recordnewStrandID

	// where &recordNewStrandID is the entire SDR data record to check if the strand is linear or circular


        // An ecDNA event touches at least two new-strand records: one linear record with the two or more flanking fragments, one
        // or more CIRCULAR records with the excised fragment(s).
        if (groupNewStrand.size() < 2)
        {
            continue;
        }

	int remainingStrandID = -1;						// newStrandID for the strand missing an excised fragment
        std::vector<SDRfragment>* flankingFragments = nullptr;			// Vector containing the pointer to the two fragments missing a central excised segment
        bool multipleRemainingRecords = false;
        std::vector<std::pair<int, std::vector<SDRfragment>*>> excisedCandidates; 					// {newStrandID, fragments}


        for (auto& [newStrandID, fragments] : groupNewStrand)
        {
            const SDRdataRecord* record = recordByNewStrand[newStrandID];	// Assign the SDR data record to a given newStrandID

            if (fragments.size() >= 2 && record->linear)			// Necessary format is one entry that is linear with two flanking fragments referencing the same oldStrandID
            {
		if (flankingFragments != nullptr)
                {
                    multipleRemainingRecords = true; // More than one candidate "remaining" record - invalid shape.
                    break;
                }

                remainingStrandID = newStrandID;
                flankingFragments = &fragments;
            }
            else if (!record->linear)			// Necessary format is one entry that is non-linear with one singular fragment referencing the same oldStrandID as the record with the two flanking fragments
            {
		// Every fragment in an ecDNA (excised) record must be acentric.
                bool hasCentromereFlag = false;

		for (const SDRfragment& fragment : fragments)
                {
                    if (fragment.hasCentromere)
                    {
                        hasCentromereFlag = true;
                        break;
                    }
                }

                if (hasCentromereFlag)
                {
                    continue; // Disqualified - not added as an excised candidate.
                }

                excisedCandidates.push_back({newStrandID, &fragments});
            }
        }


	if (multipleRemainingRecords || flankingFragments == nullptr || excisedCandidates.empty())	// If desired format not detected, skip this record
        {
            continue;
        }

	// Sort flanking fragments in ascending order of start position to find the gap between flanking fragments where the excised circular fragment was originally.
        std::sort(flankingFragments->begin(), flankingFragments->end(), [](const SDRfragment& a, const SDRfragment& b)
        {
            return a.oldStartPosition < b.oldStartPosition;
        });

	const double ecDNAtolerance = 0.001;					// tolerance in Mbp difference between contiguous segments <= 100 bp

	struct Gap
        {
            double startPos;
            double endPos;
        };

        std::vector<Gap> gaps;


	for (std::size_t i = 0; i + 1 < flankingFragments->size(); ++i)
        {
            const double gapStart = (*flankingFragments)[i].oldEndPosition;
            const double gapEnd = (*flankingFragments)[i + 1].oldStartPosition;

            if (gapEnd <= gapStart)
            {
                gaps.clear();
                break; 				// Overlap/inversion shape - invalid.
            }

            gaps.push_back({gapStart, gapEnd});
        }

	if (gaps.empty())
        {
            continue;
        }

	// Total fragments across all excised candidates must match the
        // total number of gaps - every gap must be filled by exactly one
        // fragment, and every excised fragment must fill exactly one gap.
        std::size_t totalExcisedFragments = 0;

	for (const auto& [newStrandID, fragments] : excisedCandidates)
        {
            totalExcisedFragments += fragments->size();
        }

        if (totalExcisedFragments != gaps.size())
        {
            continue;
        }

	std::vector<bool> gapClaimed(gaps.size(), false);
        std::vector<SDRecDNAevent> strandEvents;
        bool allMatched = true;

	for (const auto& [newStrandID, fragments] : excisedCandidates)
        {
            std::vector<std::pair<double, double>> matchedSegments;

            for (const SDRfragment& fragment : *fragments)
            {
                bool matched = false;

                for (std::size_t i = 0; i < gaps.size(); ++i)
                {
                    if (gapClaimed[i])
                    {
                        continue;
                    }

		    if (approxEqual(fragment.oldStartPosition, gaps[i].startPos, ecDNAtolerance) &&
                        approxEqual(fragment.oldEndPosition, gaps[i].endPos, ecDNAtolerance))
                    {
                        gapClaimed[i] = true;
                        matched = true;
                        matchedSegments.push_back({gaps[i].startPos, gaps[i].endPos});
                        break;
                    }
                }

                if (!matched)
                {
                    allMatched = false;
                    break;
                }
            }

	    if (!allMatched)
            {
                break;
            }

	    SDRecDNAevent event{};
            event.oldStrandID = oldStrandID;
            event.ecDNAsegments = matchedSegments;
            event.remainingStrandID = remainingStrandID;
            event.excisedStrandID = newStrandID;

            strandEvents.push_back(event);
        }

	if (!allMatched)
        {
            continue; // Some gap or fragment didn't match cleanly - reject the whole strand's grouping.
        }

        ecDNAevents.insert(ecDNAevents.end(), strandEvents.begin(), strandEvents.end());

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
    const SDRsubHeader& subHeader,						// Loop through the SDR subheaders to get the data records for each cell
    int numOriginalStrands)							// Determine number of original strands to determine which records are mutations
{
    std::vector<SDRdeletionInversionEvent> delInvs;

    std::map<int, std::vector<std::pair<int, SDRfragment>>> groupOldStrand;

    const double delInvTolerance = 0.001;					// ~ 0.001 Mbp = 1000 bp tolerance for position mismatch

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (!record.linear)							// If strand is non-linear, skip, not delInv
        {
            continue;
        }

        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))	// If newStrandID is not greater than the numOriginalStrands, this is not a mutation, can ignore.
        {
            continue;
        }

        for (const SDRfragment& fragment : record.fragments)
        {
            groupOldStrand[fragment.oldStrandID].push_back({record.newStrandID, fragment});	// Store the newStrandIDs and fragments making up a given oldStrandID
        }
    }

    for (auto& [oldStrandID, entries] : groupOldStrand)
    {
        std::map<int, std::vector<SDRfragment>> groupNewStrand;

        for (const auto& [newStrandID, fragment] : entries)			// Loop through the fragments associated with a given oldStrandID
        {
            groupNewStrand[newStrandID].push_back(fragment);			// Store the fragment information for each given newStrandID
        }


        if (groupNewStrand.size() != 2)						// Must touch exactly two new-strand records: one with three fragments (2 normal + 1 inverted), one with the single excised fragment.
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
            continue; 									// The excised piece itself should not be reversed (not an inversion).
        }


        int reversedCount = 0;

        for (const SDRfragment& fragment : *threeFragments)				// Check how many inversions in the remainingStrandID strand, should have only 1
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

        std::vector<NormalizedFragment> normalizedFragment;					// Fragment not necessarily structured with start and end location, but by lower position value and higher position value
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
            if (approxEqual(normalizedFragment[i].high, normalizedFragment[i + 1].low, delInvTolerance))	// Contiguous - no gap here. Skip this record
            {
                continue;
            }

            if (normalizedFragment[i + 1].low <= normalizedFragment[i].high)			// Overlapping, not a gap - invalid shape. Skip this record
            {
                overlapFound = true;
                break;
            }

            gapStart = normalizedFragment[i].high;						// If shape of mutation is valid, store fragment gap info
            gapEnd = normalizedFragment[i + 1].low;
            ++gapCount;
        }

        if (overlapFound || gapCount != 1)							// Need there to be an overlap found and exactly one gap
        {
            continue;
        }

        if (!approxEqual(excisedFragment.oldStartPosition, gapStart, delInvTolerance) ||	// If gap does not align with excised fragment within tolerance, skip this record
            !approxEqual(excisedFragment.oldEndPosition, gapEnd, delInvTolerance))
        {
            continue;
        }

        const auto reversedIt = std::find_if(normalizedFragment.begin(), normalizedFragment.end(), [](const NormalizedFragment& f)	// Look through the fragments for the inversion position
        {
            return f.reversed;
        });

        if (reversedIt == normalizedFragment.end())						// If no inversion found, skip record
        {
            continue;
        }

	// All checks passed, store deletion-inversion information
        SDRdeletionInversionEvent event{};
        event.oldStrandID = oldStrandID;
        event.inversionStart = reversedIt->low;
        event.inversionEnd = reversedIt->high;
        event.deletionStart = gapStart;
        event.deletionEnd = gapEnd;
        event.remainingStrandID = remainingStrandID;
        event.excisedStrandID = excisedStrandID;

        delInvs.push_back(event);								// Store deletion-inversion info
    }

    return delInvs;										// In writeCellDataSummary, return the size of the delInvs for the number of detected mutations
}











// Checks whether two fragments referencing the same old strand ID have a genuine gap between them (rather than meeting cleanly at a
// breakpoint) - i.e. material between them is missing. If so, populates gapStart/gapEnd with the gap's bounds.
inline bool fragmentsHaveGap(const SDRfragment& fragment1, const SDRfragment& fragment2, double& gapStart, double& gapEnd)
{
    const double delTolerance = 0.001; 				// Tolerance for base pair positioning mismatch for deletions is 0.001 Mbp = 1000 bp or less

    if (fragment1.oldStrandID != fragment2.oldStrandID)		// Both fragments being compared must have the same old strand ID
    {
        return false;
    }

    if (isReversedFragment(fragment1) || isReversedFragment(fragment2))	// Neither fragment being compared should be inverted
    {
        return false;
    }

    // If fragment 1 end is less than fragment 2 start and the positions are approximately equal, then a gap exists and compute the gap
    if (fragment1.oldEndPosition < fragment2.oldStartPosition &&
        !approxEqual(fragment1.oldEndPosition, fragment2.oldStartPosition, delTolerance))
    {
        gapStart = fragment1.oldEndPosition;
        gapEnd = fragment2.oldStartPosition;
        return true;
    }

    // If fragment 2 end is less than fragment 1 start and the positions are approximately equal, then a gap exists and compute the gap
    if (fragment2.oldEndPosition < fragment1.oldStartPosition &&
        !approxEqual(fragment2.oldEndPosition, fragment1.oldStartPosition, delTolerance))
    {
        gapStart = fragment2.oldEndPosition;
        gapEnd = fragment1.oldStartPosition;
        return true;
    }

    // If th two fragment ends are not approximately equal and one fragment's end is not less than the other fragment's start, no fragment gap
    return false;
}











// --------------------------------------------------------------- //
// Function to detect DELETION-TRANSLOCATIONS in a single cell     //
// --------------------------------------------------------------- //

// Looks for 2 fragments each, from 2 different old strands, linear, non-reversed. If BOTH strands meet
// cleanly at a breakpoint, that's a plain balanced translocation - not this. If exactly ONE strand meets cleanly and the other has a
// gap, a third record (single fragment, from the gapped strand, exactly filling the gap) confirms a deletion-translocation.

inline std::vector<SDRdeletionTranslocationEvent> detectDeletionTranslocations(
    const SDRsubHeader& subHeader,						// Read SDR subheader to loop through each cell's damage record
    int numOriginalStrands)							// numOriginalStrands tells us which data records are mutations
{
    std::vector<SDRdeletionTranslocationEvent> delTras;

    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;

    const double delTraTolerance = 0.001;					// Tolerance for the mismathc in contiguous base positions is 0.001 Mbp = 1000 bp or less

    for (std::size_t i = 0; i < records.size(); ++i)
    {
        SDRfragment fragmentA1{};
        SDRfragment fragmentB1{};

	// If fragments A1 and B1 are not possible translocations, skip this data record
        if (!getTranslocationCandidateFragments(records[i], numOriginalStrands, fragmentA1, fragmentB1))
        {
            continue;
        }

	// Checking the proceeding record in the data section if the translocation is present there instead
        for (std::size_t j = i + 1; j < records.size(); ++j)
        {
            SDRfragment fragmentA2{};
            SDRfragment fragmentB2{};

	    // If fragments A2 and B2 are not possible translocations, skip this data record
            if (!getTranslocationCandidateFragments(records[j], numOriginalStrands, fragmentA2, fragmentB2))
            {
                continue;
            }

            SDRfragment firstMatch{};
            SDRfragment secondMatch{};


	    // Find which fragments match up and distinguish from the deleted portion
	    // This section determines which fragments should be compared in fragmentsShareBreakpoint and fragmentsHaveGap functions
	    // Comparing fragment old strand IDs across two consecutive data records (the first record should have fragments A1 and B1, 
	    // then the second record should have fragments A2 and B2)
            if (fragmentA2.oldStrandID == fragmentA1.oldStrandID && fragmentB2.oldStrandID == fragmentB1.oldStrandID)
            {
                firstMatch = fragmentA2;
                secondMatch = fragmentB2;
            }
            else if (fragmentA2.oldStrandID == fragmentB1.oldStrandID && fragmentB2.oldStrandID == fragmentA1.oldStrandID)
            {
                firstMatch = fragmentB2;
                secondMatch = fragmentA2;
            }
            else
            {
                continue;
            }

	    // Determine the breakpoint locations for the two translocated records
            double breakpointA = 0.0;
            const bool cleanBreakA = fragmentsShareBreakpoint(fragmentA1, firstMatch, breakpointA);

            double breakpointB = 0.0;
            const bool cleanBreakB = fragmentsShareBreakpoint(fragmentB1, secondMatch, breakpointB);

            if (cleanBreakA && cleanBreakB)
            {
                continue; 				// Plain balanced translocation - handled elsewhere, not here.
            }

	    // Determine the deleted section length since not a clean breakpoint (not balanced translocation)
            double gapStart = 0.0;					// Start position of the deleted segment in Mbp
            double gapEnd = 0.0;					// End position of the deleted segment in Mbp
            int deletedOldStrandID = -1;				// Track origin of the deleted segment
            double cleanBreakPos = 0.0;					// The deleted segment should have a clean breakpoint with one of the translocated fragment ends

	    // Check which original strand contains the deleted segment
	    // If the deletion occurs in cleanBreakB, then the fragments B1 and secondMatch should have a gap
            if (cleanBreakA && !cleanBreakB)
            {
                if (!fragmentsHaveGap(fragmentB1, secondMatch, gapStart, gapEnd))	// If there's no gap, there is no deletion, skip to next record
                {
                    continue;
                }

                deletedOldStrandID = fragmentB1.oldStrandID;				// Store the original strand ID of the deleted segment
                cleanBreakPos = breakpointA;						// Store the location of the break of the deletion
            }
	    // If the deletion occurs in cleanBreakA, then the fragments A1 and firstMatch should have a gap
            else if (cleanBreakB && !cleanBreakA)
            {
                if (!fragmentsHaveGap(fragmentA1, firstMatch, gapStart, gapEnd))
                {
                    continue;
                }

                deletedOldStrandID = fragmentA1.oldStrandID;				// Store the original strand ID of the deleted segment
                cleanBreakPos = breakpointB;						// Store the location of the break of the deletion
            }
            else
            {
                continue; 								// Neither strand cleanly matches - not this mutation type.
            }

	    // Now looking at the deletion record of the deletion-translocation
            for (const SDRdataRecord& deletionRecord : records)
            {
		// Deleted portion should have its own unique newStrandID, because it is a new data entry, otherwise, skip this record
                if (deletionRecord.newStrandID == records[i].newStrandID || deletionRecord.newStrandID == records[j].newStrandID)
                {
                    continue;
                }

		// Deleted portion must be linear, otherwise not a deletion-translocation
                if (!deletionRecord.linear)
                {
                    continue;
                }

		// Safety check, make sure the deletion containing a single fragment entry has a newStrandID > numOriginalStrands, otherwise it is an
		// intact strand entry
                if (!isRearrangementCandidate(deletionRecord.newStrandID, numOriginalStrands))
                {
                    continue;
                }

		// Deletion record must contain one fragment
                if (deletionRecord.fragments.size() != 1)
                {
                    continue;
                }

		// Confirmed the excised fragment is the deletion
                const SDRfragment& excisedFragment = deletionRecord.fragments[0];

		// Deleted fragment's old strand ID should correspond to one of the other records' old strand IDs that underwent the translocation
                if (excisedFragment.oldStrandID != deletedOldStrandID)
                {
                    continue;
                }

		// The deleted section is not inverted, otherwise skip this record
                if (isReversedFragment(excisedFragment))
                {
                    continue;
                }

		// The start and end locations of the excised fragment should be approximately equal to the gap left on the original strand that underwent
		// a translocation (within 1000 bp tolerance)
                if (!approxEqual(excisedFragment.oldStartPosition, gapStart, delTraTolerance) ||
                    !approxEqual(excisedFragment.oldEndPosition, gapEnd, delTraTolerance))
                {
                    continue;
                }

		// All checks for deletion-translocation passed, store information in event object
                SDRdeletionTranslocationEvent event{};
                event.oldStrandA = fragmentA1.oldStrandID;
                event.oldStrandB = fragmentB1.oldStrandID;
                event.deletedOldStrandID = deletedOldStrandID;
                event.cleanBreakPos = cleanBreakPos;
                event.deletionStart = gapStart;
                event.deletionEnd = gapEnd;
                event.newStrandID1 = records[i].newStrandID;
                event.newStrandID2 = records[j].newStrandID;
                event.excisedStrandID = deletionRecord.newStrandID;

                delTras.push_back(event);						// Store deletion-translocation information
                break;
            }
        }
    }

    return delTras;									// In writeCellDataSummary, return the size of the delTras vector to summarize the number of mutations
}











#endif
