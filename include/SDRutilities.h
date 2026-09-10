#ifndef SDR_UTILITIES_H
#define SDR_UTILITIES_H

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>
#include <set>

#include "SDRtypes.h"

// -------------------------------------------------- //
// --------- SDR MUTATION DETECTION HELPERS --------- //
// -------------------------------------------------- //

// oldStrandID always refers to one of the
// original intact chromosomes (46 for humans), never to a
// newly-created strand ID produced by an earlier rearrangement.
// This is what allows every detect<Mutation>() function below to group a
// cell's fragments by oldStrandID and treat that as "all the
// rearrangement activity affecting a given original chromosome",
// rather than needing to trace strand lineage across events.



// Determines which original strand a record "belongs to" for karyogram
// slot placement. A record with exactly one centromere-bearing
// fragment unambiguously belongs to that fragment's strand - this is
// the correct rule for any BALANCED rearrangement, since by
// definition each derivative keeps exactly one parent's centromere.
//
// FUTURE WORK: a record with zero (acentric) or two+ (dicentric)
// centromere-bearing fragments currently falls back to the old
// first-fragment rule, placing it in exactly one slot. True dicentric
// support would need this record to be placed in MULTIPLE slots at
// once (one per centromere it carries) - that requires changing
// recordsByOriginalStrand's grouping itself (a record can currently
// only live under one key), not just this function. Left as a single
// deliberate fallback here so that future change has one clear call
// site to revisit, rather than being scattered.
inline int determineHomeStrandID(const SDRdataRecord& record)
{
    int centromereStrandID = -1;
    int centromereCount = 0;

    for (const SDRfragment& fragment : record.fragments)
    {
        if (fragment.hasCentromere)
        {
            centromereStrandID = fragment.oldStrandID;
            ++centromereCount;
        }
    }

    if (centromereCount == 1)
    {
        return centromereStrandID;
    }

    // Zero or 2+ centromeres - ambiguous under today's single-slot
    // model. Fall back to the original strand-ordering convention.
    return record.fragments.empty() ? -1 : record.fragments[0].oldStrandID;
}









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


    // Tracks each newStrandID's TRUE total fragment count across its whole
    // record - used to reject a candidate whose fragments are only
    // PARTLY from the old strand being examined (e.g. a deletion-insertion
    // recipient record, which mixes native fragments from one old strand
    // with one foreign fragment from another).
    std::map<int, std::size_t> newStrandIDtotalFragments;


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

        newStrandIDtotalFragments[record.newStrandID] = record.fragments.size();


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
	    if (fragments.size() != newStrandIDtotalFragments[newStrandID]) // Reject if this record has fragments from OTHER old strands too - a pure remaining/excised piece must be entirely from this old strand.
            {
                continue;
            }

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
	// Generalized check for strands with multiple deletions and gaps.
        std::vector<bool> excisedClaimed(excisedCandidates.size(), false);	// Store whether the gaps in a record were found by fragments in other records
        std::vector<SDRdeletionEvent> deletionEvents;				// Store the relevant information regarding strand IDs and fragments sizes
        bool allGapsMatched = true;

	for (const Gap& gap : gaps)						// Loop through all the gaps in a given data record (two flanking fragments that are not contiguous within tolerance of ~1000 bp)
        {
            bool matched = false;

            for (std::size_t i = 0; i < excisedCandidates.size(); ++i)
            {
                if (excisedClaimed[i])						// If the excised fragment has not be claimed (false) skip this record
                {
                    continue;
                }

                const SDRfragment& excisedFragment = *excisedCandidates[i].second;	// Access second member of pair (i.e. the specific fragment)

                // Use a tolerance of ~1000 bases for a long deletion event to be tracked.
		// If the gap ends match up roughly with the fragment ends, then claim this excised fragment and the gap and fragment match
                if (approxEqual(excisedFragment.oldStartPosition, gap.startPos, delTolerance) &&
                    approxEqual(excisedFragment.oldEndPosition, gap.endPos, delTolerance))
                {
                    excisedClaimed[i] = true;
                    matched = true;

                    // Store the locations of the deletion and the strands involved for later Karyogram plotting
                    SDRdeletionEvent event{};					// Construct the SDRdeletionEvent object 'event'
		    event.oldStrandID = oldStrandID;				// Store this record's old strand ID (all fragments should have the same old strand ID)
                    event.deletionStart = gap.startPos;				// Store the start position of the deletion in Mbp
                    event.deletionEnd = gap.endPos;				// Store the end position of the deletion in Mbp
                    event.remainingStrandID = remainingStrandID;		// Store the new strand ID of the record containing gaps
                    event.excisedStrandID = excisedCandidates[i].first;		// Store the new strand ID of the record containing the excised fragment

                    deletionEvents.push_back(event);				// Store the SDRdeletionEvent object in the deletionEvents vector
                    break;
                }
            }

            if (!matched)
            {
                allGapsMatched = false;						// All gaps must match to have multiple long deletions in a single original strand, otherwise this is some other more complex mutation type.
                break;
            }
        }

	if (!allGapsMatched)                                                    // Some gap didn't match any excised record - not a valid deletion shape.
        {
            continue;
        }

        deletions.insert(deletions.end(), deletionEvents.begin(), deletionEvents.end());	// Append to the deletions vector if all checks have been satisfied.
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
// intact region with no gaps (balanced) and can have more than one reverse fragment.
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

            if (reversedCount < 1)						// Need at least one reversal for a balanced inversion to be detected
            {
                continue;
            }

            struct NormalizedFragment						// Normalize fragment position to account for potential old strand start positions being greater than old strand end positions
            {
                double lowerPos;							// The lower fragment start position (closer to the start of the p arm)
                double higherPos;							// The higher fragment end position (closer to the end of the q arm)
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
                    return a.lowerPos < b.lowerPos;
                });

	    // Contiguous if the ends of fragments match within ~10 bp
            bool contiguous = true;						// Minimal loss of bases (<= 10) can occur in balanced inversion, contiguous must be true
            for (std::size_t i = 0; i + 1 < normalized.size(); i++)		// Loop through all normalized fragments per old Strand ID
            {
                if (!approxEqual(normalized[i].higherPos, normalized[i + 1].lowerPos, balInvTolerance))
                {
                    contiguous = false;
                    break;
                }
            }

            if (!contiguous)
            {
                continue;
            }


	    // Emit one event per reversed fragment found - N reversed
            // fragments means N separate inversions on this strand.
            for (const NormalizedFragment& fragment : normalized)
            {
                if (!fragment.reversed)
                {
                    continue;
                }

                // Record inversion event details for karyogram plotting later
                SDRinversionEvent event{};
                event.oldStrandID = oldStrandID;
                event.inversionStart = fragment.lowerPos;
                event.inversionEnd = fragment.higherPos;
                event.newStrandID = record.newStrandID;

                inversions.push_back(event);                          		     	// Append inversion event details to the inversions vector
            }
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


// Generalized to detect translocations anywhere on a chromosome (not just at its ends),
// and to correctly handle an exchanged piece that re-attaches in reversed orientation
// (which looks locally like an inversion but isn't one, since it spans two strands).
// Pools every fragment referencing a given old strand across BOTH candidate records,
// normalizes each to [min, max] (so reversed fragments sort correctly).
inline std::vector<SDRtranslocationEvent> detectTranslocations(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader)
{
    std::vector<SDRtranslocationEvent> translocations;					// Translocations vector to store translocation events

    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;			// Records vector to store the SDR data records/entries

    const double balTraTolerance = 0.00001;


    auto getRecordStrandFragments = [&](const SDRdataRecord& record, std::map<int, std::vector<SDRfragment>>& groupStrand) -> bool
    {
        if (!record.linear)
        {
            return false;
        }

        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))
        {
            return false;
        }

        for (const SDRfragment& fragment : record.fragments)
        {
            groupStrand[fragment.oldStrandID].push_back(fragment);
        }

        return !groupStrand.empty();
    };



    // Checks that the pooled fragments for one strand, normalized to
    // [min,max] and sorted, tile [0, fullSize] exactly with no gaps or
    // overlaps. On success, breakpointsOut holds every internal
    // boundary between consecutive fragments.
    auto checkFullTiling = [&](std::vector<SDRfragment> pooled, int strandID, std::vector<double>& breakpointsOut) -> bool
    {
        const std::size_t sizeIndex = static_cast<std::size_t>(strandID);

        if (sizeIndex >= masterHeader.intactChromosomeSizes.size())
        {
            return false;
        }

        const double fullSize = masterHeader.intactChromosomeSizes[sizeIndex];

        if (fullSize <= 0.0)
        {
            return false;
	}


	struct NormalizedFragment
        {
            double low;
            double high;
        };

        std::vector<NormalizedFragment> normalized;
        normalized.reserve(pooled.size());

	for (const SDRfragment& fragment : pooled)
        {
            normalized.push_back({
                std::min(fragment.oldStartPosition, fragment.oldEndPosition),
                std::max(fragment.oldStartPosition, fragment.oldEndPosition)
            });
        }

        std::sort(normalized.begin(), normalized.end(), [](const NormalizedFragment& a, const NormalizedFragment& b)
        {
            return a.low < b.low;
        });

        if (!approxEqual(normalized.front().low, 0.0, balTraTolerance))
        {
            return false;
        }

	if (!approxEqual(normalized.back().high, fullSize, balTraTolerance))
        {
            return false;
        }

        breakpointsOut.clear();

        for (std::size_t k = 0; k + 1 < normalized.size(); ++k)
        {
            if (!approxEqual(normalized[k].high, normalized[k + 1].low, balTraTolerance))
            {
                return false;
            }

            breakpointsOut.push_back(normalized[k].high);
        }

        return true;
    };


    for (std::size_t i = 0; i < records.size(); ++i)
    {
        std::map<int, std::vector<SDRfragment>> iByStrand;

        if (!getRecordStrandFragments(records[i], iByStrand))
        {
            continue;
        }

        if (iByStrand.size() != 2)
        {
            continue;
        }

	for (std::size_t j = i + 1; j < records.size(); ++j)
        {
            std::map<int, std::vector<SDRfragment>> jByStrand;

            if (!getRecordStrandFragments(records[j], jByStrand))
            {
                continue;
            }

            if (jByStrand.size() != 2)
            {
                continue;
            }

            std::vector<int> iStrands;


	    for (auto& [strandID, fragments] : iByStrand)
            {
                iStrands.push_back(strandID);
            }

            std::vector<int> jStrands;

            for (auto& [strandID, fragments] : jByStrand)
            {
                jStrands.push_back(strandID);
            }

            if (iStrands != jStrands)
            {
                continue;
            }

	    const int strandA = iStrands[0];
            const int strandB = iStrands[1];

            std::vector<SDRfragment> pooledA = iByStrand[strandA];
            pooledA.insert(pooledA.end(), jByStrand[strandA].begin(), jByStrand[strandA].end());

            std::vector<SDRfragment> pooledB = iByStrand[strandB];
            pooledB.insert(pooledB.end(), jByStrand[strandB].begin(), jByStrand[strandB].end());

            std::vector<double> breakpointsA;
            std::vector<double> breakpointsB;

            if (!checkFullTiling(pooledA, strandA, breakpointsA) || !checkFullTiling(pooledB, strandB, breakpointsB))
            {
                continue;
            }


	    SDRtranslocationEvent event{};
            event.oldStrandA = strandA;
            event.oldStrandB = strandB;
            event.breakpointsA = breakpointsA;
            event.breakpointsB = breakpointsB;
            event.newStrandID1 = records[i].newStrandID;
            event.newStrandID2 = records[j].newStrandID;

            translocations.push_back(event);
	}
    }

    return translocations;
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


    // Tracks each newStrandID's TRUE total fragment count across its whole
    // record - used below to reject a candidate whose fragments are only
    // PARTLY from the old strand being examined.
    std::map<int, std::size_t> newStrandIDtotalFragments;


    for (const SDRdataRecord& record : subHeader.dataRecords)			// Loop through each data record per cell subheader
    {
	if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))	// Check if a data record consists of a rearranged fragment, if not, skip
	{
	    continue;
	}

        newStrandIDtotalFragments[record.newStrandID] = record.fragments.size();


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

	    if (fragments.size() != newStrandIDtotalFragments[newStrandID]) // Reject if this record has fragments from OTHER old strands too.
            {
                continue;
            }

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

	struct Gap								// Gap struct to store the start and end positions of multiple gaps caused by multiple ecDNA events in a single fragment
        {
            double startPos;
            double endPos;
        };

        std::vector<Gap> gaps;


	for (std::size_t i = 0; i + 1 < flankingFragments->size(); ++i)		// Loop through all flanking fragments, should be a gap between each
        {
            const double gapStart = (*flankingFragments)[i].oldEndPosition;	// Gap start position in Mbp
            const double gapEnd = (*flankingFragments)[i + 1].oldStartPosition;	// Gap end position in Mbp

            if (gapEnd <= gapStart)						// Overlap/inversion shape - invalid.
            {
                gaps.clear();
                break;
            }

            gaps.push_back({gapStart, gapEnd});					// If valid shape, store the gap positions in the gap vector
        }

	if (gaps.empty())							// If no valid gaps found, skip this record entry.
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

        if (totalExcisedFragments != gaps.size())				// Number of excised fragments must equal number of gaps for multiple ecDNA
        {
            continue;
        }

	std::vector<bool> gapClaimed(gaps.size(), false);			// Initialize all the vector as false by default
        std::vector<SDRecDNAevent> strandEvents;				// Append to strandEvents if gaps match excised fragments and flanking fragments
        bool allMatched = true;							// All gaps must match to have multiple ecDNA on the same strand

	for (const auto& [newStrandID, fragments] : excisedCandidates)		// Check all excised candidates
        {
            std::vector<std::pair<double, double>> matchedSegments;		// Store the start and end positions of the gaps

            for (const SDRfragment& fragment : *fragments)			// Check the flanking fragment start and end positions and compare with excised fragments
            {
                bool matched = false;

                for (std::size_t i = 0; i < gaps.size(); ++i)			// Loop through the number of gaps in a given data record
                {
                    if (gapClaimed[i])						// If a gap was not found, skip this record
                    {
                        continue;
                    }

		    // If the ends of the flanking fragments match the ends of the gaps within tolerance, then this gap has been claimed and matches
		    // An excised fragment
		    if (approxEqual(fragment.oldStartPosition, gaps[i].startPos, ecDNAtolerance) &&
                        approxEqual(fragment.oldEndPosition, gaps[i].endPos, ecDNAtolerance))
                    {
                        gapClaimed[i] = true;
                        matched = true;
                        matchedSegments.push_back({gaps[i].startPos, gaps[i].endPos});
                        break;
                    }
                }

                if (!matched)							// All gaps must match excised fragments, otherwise this is a more complex mutation, not multiple ecDNA in a single strand.
                {
                    allMatched = false;
                    break;
                }
            }

	    if (!allMatched)
            {
                break;
            }

	    // All checks for mutation shape passed, store the mutation information.
	    SDRecDNAevent event{};
            event.oldStrandID = oldStrandID;
            event.ecDNAsegments = matchedSegments;
            event.remainingStrandID = remainingStrandID;
            event.excisedStrandID = newStrandID;

            strandEvents.push_back(event);
        }

	if (!allMatched)							// Some gap or fragment didn't match cleanly - reject the whole strand's grouping.
        {
            continue;
        }


	// Append this ecDNA event to the ecDNAevents vector
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
















// ------------------------------------------------------------------------------------ //
// Function to detect DELETION-INSERTION entries in SDR file
// ------------------------------------------------------------------------------------ //

// A deletion-insertion event is characterized by two new-strand records:
// 1) A "donor" record: two fragments from the same old strand ID, with a real gap
// between them (identical shape to a plain long deletion's "remaining" piece).
// 2) A "recipient" record: three fragments - two from a DIFFERENT old strand ID that
// are themselves contiguous (the recipient's own material, split by the insertion),
// and one FOREIGN fragment (from the donor's old strand ID) sitting between them,
// whose position exactly matches the donor's gap.
// Unlike a balanced translocation, only ONE segment moves - the donor loses material
// with nothing coming back, and the recipient's own material is fully retained, just
// split by the inserted piece.
inline std::vector<SDRdeletionInsertionEvent> detectDeletionInsertions(const SDRsubHeader& subHeader, int numOriginalStrands)
{
    double tolerance = 0.001;                                              	// ~1000 bp tolerance, matching the long deletion convention.

    std::vector<SDRdeletionInsertionEvent> delInsEvents;			// Store the deletion insertion mutation strand and fragment information

    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;		// SDR data records of a given cell

    for (std::size_t i = 0; i < records.size(); i++)                       	// Outer loop: look for candidate DONOR records (strands with a deletion)
    {
        const SDRdataRecord& donorRecord = records[i];

        if (!donorRecord.linear)						// Cannot have circular fragments in deletion-insertions
        {
            continue;
        }

        if (!isRearrangementCandidate(donorRecord.newStrandID, numOriginalStrands))	// Check if this record contains a mutation, if not then skip
        {
            continue;
        }

        if (donorRecord.fragments.size() != 2)                             	// Donor "remaining" piece must have exactly two fragments
        {
            continue;
        }

        const SDRfragment& donorFragA = donorRecord.fragments[0];		// Donor will have two remaining fragments flanking the deleted segment
        const SDRfragment& donorFragB = donorRecord.fragments[1];

        if (donorFragA.oldStrandID != donorFragB.oldStrandID)			// Both flanking fragments must be from the same original strand
        {
            continue;
        }

        if (isReversedFragment(donorFragA) || isReversedFragment(donorFragB))	// Neither flanking fragment can be inverted, otherwise this is a more complex mutation
        {
            continue;
        }

        const int donorOldStrandID = donorFragA.oldStrandID;			// Same as donorFragB.oldStrandID

        SDRfragment lowerDonorFrag = donorFragA;				// Flanking fragment on donor with the lower start position
        SDRfragment upperDonorFrag = donorFragB;				// Flanking fragment on donor with the higher start position

        if (lowerDonorFrag.oldStartPosition > upperDonorFrag.oldStartPosition)
        {
            std::swap(lowerDonorFrag, upperDonorFrag);
        }

	// Gap start and end positions must fit between the flaning fragment ends
        const double segmentStart = lowerDonorFrag.oldEndPosition;
        const double segmentEnd = upperDonorFrag.oldStartPosition;

        if (segmentEnd <= segmentStart)                                    	// Must be an actual gap
        {
            continue;
        }

        for (std::size_t j = 0; j < records.size(); j++)                   	// Inner loop: look for a matching candidate RECIPIENT record
        {
            if (i == j)
            {
                continue;
            }

            const SDRdataRecord& recipientRecord = records[j];

            if (!recipientRecord.linear)					// All strands must be linear in deletion-insertions
            {
                continue;
            }

            if (!isRearrangementCandidate(recipientRecord.newStrandID, numOriginalStrands))	// Make sure the record contains a mutation, otherwise skip
            {
                continue;
            }

            if (recipientRecord.fragments.size() != 3)                     	// Recipient: two native + one foreign insert
            {
                continue;
            }

            bool anyReversed = false;						// Check for inversions, should not have any

            for (const SDRfragment& fragment : recipientRecord.fragments)
            {
                if (isReversedFragment(fragment))				// Cannot have reverse fragments, this is a more complex mutation
                {
                    anyReversed = true;
                    break;
                }
            }

            if (anyReversed)
            {
                continue;
            }

            std::vector<SDRfragment> nativeFragments;				// Native fragments stay on the original strand, do not transfer to another strand
            const SDRfragment* foreignFragment = nullptr;
            int recipientOldStrandID = -1;

            for (const SDRfragment& fragment : recipientRecord.fragments)	// Check all recipient records (3 fragment records with 2 of the same old strand IDs and 1 foreign old strand ID)
            {
                if (fragment.oldStrandID == donorOldStrandID)              	// Foreign fragment must trace back to the donor's old strand
                {
                    if (foreignFragment != nullptr)                        	// Only one foreign fragment allowed
                    {
                        foreignFragment = nullptr;
                        break;
                    }

                    foreignFragment = &fragment;				// Foreign fragment found
                }
                else
                {
                    if (recipientOldStrandID == -1)
                    {
                        recipientOldStrandID = fragment.oldStrandID;
                    }
                    else if (fragment.oldStrandID != recipientOldStrandID) 	// Native fragments must all share the same old strand ID
                    {
                        recipientOldStrandID = -1;
                        break;
                    }

                    nativeFragments.push_back(fragment);			// Native fragment found
                }
            }

            if (foreignFragment == nullptr || recipientOldStrandID == -1 || nativeFragments.size() != 2)
            {
                continue;
            }

	    // The gap left by the deletion must match within tolerance the ends of the flanking fragments
            if (!approxEqual(foreignFragment->oldStartPosition, segmentStart, tolerance) ||
                !approxEqual(foreignFragment->oldEndPosition, segmentEnd, tolerance))
            {
                continue;                                                  	// Foreign fragment must match the donor's deleted segment exactly
            }

            SDRfragment lowerNativeFrag = nativeFragments[0];			// Native fragment with the lower start position
            SDRfragment upperNativeFrag = nativeFragments[1];			// Native fragment with the higher start position

            if (lowerNativeFrag.oldStartPosition > upperNativeFrag.oldStartPosition)
            {
                std::swap(lowerNativeFrag, upperNativeFrag);
            }

	    // Recipient's own material must be cleanly split, nothing lost
            if (!approxEqual(lowerNativeFrag.oldEndPosition, upperNativeFrag.oldStartPosition, tolerance))
            {
                continue;
            }

	    // Insertion occurs at the end of the first flanking fragment
            const double insertionPoint = lowerNativeFrag.oldEndPosition;

	    // All checks passed for deletion-insertion shape, store appropriate mutation information
            SDRdeletionInsertionEvent event{};
            event.donorOldStrandID = donorOldStrandID;
            event.recipientOldStrandID = recipientOldStrandID;
            event.segmentStart = segmentStart;
            event.segmentEnd = segmentEnd;
            event.insertionPoint = insertionPoint;
            event.donorRemainingNewStrandID = donorRecord.newStrandID;
            event.recipientNewStrandID = recipientRecord.newStrandID;

            delInsEvents.push_back(event);
        }
    }

    return delInsEvents;
}











// ------------------------------------------------------------------------------------ //
// Function to detect CHROMOPLEXY events in SDR file
// ------------------------------------------------------------------------------------ //

// Builds a graph where each old strand ID is a node. For every rearranged record, every
// PAIR of distinct old strand IDs referenced by its fragments gets an edge - regardless
// of whether that record matches any specific named mutation shape (translocation,
// deletion-translocation, etc.). Any connected component spanning 3 or more distinct
// strands is reported as one chromoplexy event.
inline std::vector<SDRchromoplexyEvent> detectChromoplexy(const SDRsubHeader& subHeader, int numOriginalStrands)
{
    std::vector<SDRchromoplexyEvent> chromoplexyEvents;

    const std::vector<SDRdeletionEvent> deletions = detectDeletions(subHeader, numOriginalStrands);


    // Initial State: Every chromosome is isolated in its own group (parent[x] = x).
    // Mutation 1 (Chromosome 1 and 2 swap material): You call unite(1, 2). Chromosome 1 and 2 are now linked. parent[1] becomes 2.
    // Mutation 2 (Chromosome 2 and 3 swap material): You call unite(2, 3). Now, Chromosome 1, 2, and 3 are all part of the same group led by 3.
    // If you later ask findRoot(1), the code will climb up (1 -> 2 -> 3), optimize the path, and return 3.

    // Groups chromosomes based on their shared mutations
    std::map<int, int> parent;

    auto findRoot = [&](int x) -> int		// x are the individual chromosomes
    {
        if (!parent.count(x))			// If chromosome x was found in the parent collection or not, if not initialize x as the parent
        {
            parent[x] = x;
        }

        while (parent[x] != x)			// Find root element
        {
            parent[x] = parent[parent[x]];	// path compression
            x = parent[x];			// Move pointer to newly assigned parent
        }

        return x;				// Return the root representative of the set of chromosomes
    };


    // Merge two chromosome mutation groups into a single combined group
    auto unite = [&](int a, int b)
    {
        const int rootA = findRoot(a);
        const int rootB = findRoot(b);

        if (rootA != rootB)			// Check if a and b belong to different groups
        {
            parent[rootA] = rootB;		// Unite the separate groups if a and b are not joined already
        }
    };


    // For every rearranged record, track which distinct old strand IDs
    // it references - and union every pair found together.
    std::vector<std::pair<int, std::set<int>>> recordStrandSets; // {newStrandID, distinct old strand IDs}

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))
        {
            continue;
        }

	std::set<int> strandsInRecord;

        for (const SDRfragment& fragment : record.fragments)
        {
            strandsInRecord.insert(fragment.oldStrandID);
        }

        if (strandsInRecord.size() < 2)
        {
            continue; // Purely intra-strand record - contributes no edges.
        }

        recordStrandSets.push_back({record.newStrandID, strandsInRecord});

        std::vector<int> strandVec(strandsInRecord.begin(), strandsInRecord.end());

	for (std::size_t a = 0; a < strandVec.size(); ++a)
        {
            for (std::size_t b = a + 1; b < strandVec.size(); ++b)
            {
                unite(strandVec[a], strandVec[b]);
            }
        }
    }

    std::map<int, std::set<int>> componentsByRoot;

    for (auto& [strandID, strandParent] : parent)
    {
        componentsByRoot[findRoot(strandID)].insert(strandID);
    }

    for (auto& [root, strandSet] : componentsByRoot)
    {
        if (strandSet.size() < 3)
        {
            continue;
        }

        std::vector<int> contributingIDs;

        for (auto& [newStrandID, strands] : recordStrandSets)
        {
            for (int s : strands)
            {
                if (strandSet.count(s))
                {
                    contributingIDs.push_back(newStrandID);
                    break;
                }
            }
        }


        SDRchromoplexyEvent event{};
        event.involvedStrandIDs = std::vector<int>(strandSet.begin(), strandSet.end());
        event.contributingNewStrandIDs = contributingIDs;

        for (const SDRdeletionEvent& deletionEvent : deletions)
        {
            if (strandSet.count(deletionEvent.oldStrandID))
            {
                event.strandsWithAccompanyingDeletions.push_back(deletionEvent.oldStrandID);
            }
        }

        chromoplexyEvents.push_back(event);
    }

    return chromoplexyEvents;

}













// ------------------------------------------------------------------------------------ //
// Function to detect CHROMOTHRIPSIS events in SDR file
// ------------------------------------------------------------------------------------ //

// Same underlying strand-connectivity graph as detectChromoplexy(), built from raw
// fragment data, so clustering doesn't depend on any named mutation shape. Components of
// size 3+ are skipped (chromoplexy's territory). For components of size 1 or 2, rather
// than counting named detected mutation events , this counts the TOTAL number of data-field-3 fragments across every
// rearranged record touching the cluster's strand(s) - a direct measure of how many
// pieces the strand(s) were actually broken into, regardless of shape. Only counts the 
// remaining fragments in the original strand, not the excised pieces, as a
// clearer metric for the number of rearrangements that occurred in a strand.
inline std::vector<SDRchromothripsisEvent> detectChromothripsis(const SDRsubHeader& subHeader, int numOriginalStrands, int minFragmentCount = 10)
{
    std::vector<SDRchromothripsisEvent> chromothripsisEvents;

    std::map<int, int> parent;

    auto findRoot = [&](int x) -> int
    {
        if (!parent.count(x))
        {
            parent[x] = x;
        }

        while (parent[x] != x)
        {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }

        return x;
    };

    auto unite = [&](int a, int b)
    {
        const int rootA = findRoot(a);
        const int rootB = findRoot(b);

        if (rootA != rootB)
        {
            parent[rootA] = rootB;
        }
    };

    // Map from old strand ID -> set of newStrandIDs (records) whose
    // fragments touch it, whether intra-strand-only or inter-strand.
    std::map<int, std::set<int>> recordsByStrand;

    // Map from newStrandID -> the record's own total fragment count, so
    // each record's fragments are only ever tallied once even if it
    // touches multiple strands in a cluster.
    std::map<int, std::size_t> fragmentCountByNewStrandID;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))
        {
            continue;
        }

        std::set<int> strandsInRecord;

        for (const SDRfragment& fragment : record.fragments)
        {
            strandsInRecord.insert(fragment.oldStrandID);
        }

        fragmentCountByNewStrandID[record.newStrandID] = record.fragments.size();

        for (int strandID : strandsInRecord)
        {
            findRoot(strandID); // Register the node even if this record touches only one strand.
            recordsByStrand[strandID].insert(record.newStrandID);
        }

        if (strandsInRecord.size() < 2)
        {
            continue; // Intra-strand record - no edge to add.
        }

        std::vector<int> strandVec(strandsInRecord.begin(), strandsInRecord.end());

        for (std::size_t a = 0; a < strandVec.size(); ++a)
        {
            for (std::size_t b = a + 1; b < strandVec.size(); ++b)
            {
                unite(strandVec[a], strandVec[b]);
            }
        }
    }

    std::map<int, std::set<int>> componentsByRoot;

    for (auto& [strandID, strandParent] : parent)
    {
        componentsByRoot[findRoot(strandID)].insert(strandID);
    }


    for (auto& [root, strandSet] : componentsByRoot)
    {
        if (strandSet.size() > 2)
        {
            continue; // 3+ strands - chromoplexy's territory, not chromothripsis.
        }

        // Any record with 2+ fragments is a "major" surviving piece and
        // counts in full - whether it's a single-strand remaining piece
        // (like a deletion's) or a recombined multi-strand derivative
        // (like a translocation/chromothripsis product). A record with
        // exactly 1 fragment is always a trivial excised byproduct
        // (true for every excised/circular shape our detectors define)
        // and gets excluded.
        std::set<int> majorNewStrandIDs;

        for (int strandID : strandSet)
        {
            for (int newStrandID : recordsByStrand[strandID])
            {
                if (fragmentCountByNewStrandID[newStrandID] >= 2)
                {
                    majorNewStrandIDs.insert(newStrandID);
                }
            }
        }

        int totalFragmentCount = 0;

        for (int newStrandID : majorNewStrandIDs)
        {
            totalFragmentCount += static_cast<int>(fragmentCountByNewStrandID[newStrandID]);
        }

        if (totalFragmentCount < minFragmentCount)
        {
            continue;
        }

        SDRchromothripsisEvent event{};
        event.involvedStrandIDs = std::vector<int>(strandSet.begin(), strandSet.end());
        event.contributingNewStrandIDs = std::vector<int>(majorNewStrandIDs.begin(), majorNewStrandIDs.end());
        event.totalFragmentCount = totalFragmentCount;

        chromothripsisEvents.push_back(event);
    }

    return chromothripsisEvents;
}










#endif
