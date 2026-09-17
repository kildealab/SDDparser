#ifndef SDR_UTILITIES_H
#define SDR_UTILITIES_H

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>
#include <set>

#include "SDRtypes.h"


namespace sddparser
{

// -------------------------------------------------- //
// --------- SDR MUTATION DETECTION HELPERS --------- //
// -------------------------------------------------- //

// oldStrandID always refers to one of the  original intact chromosomes (46 for humans), never to a
// newly-created strand ID produced by an earlier rearrangement. This is what allows every detect<Mutation>() 
// function below to group a cell's fragments by oldStrandID and treat that as "all the rearrangement activity 
// affecting a given original chromosome", rather than needing to trace strand lineage across events.


// Determines which original strand a record "belongs to" for karyogram slot placement. A record with exactly one centromere-bearing
// fragment unambiguously belongs to that fragment's strand, this is the correct rule for any BALANCED rearrangement, since by
// definition each derivative keeps exactly one parent's centromere.
inline int determineHomeStrandID(const SDRdataRecord& record)
{
    int centromereStrandID = -1;
    int centromereCount = 0;

    for (const SDRfragment& fragment : record.fragments)			// Find which fragment holds the centromere
    {
        if (fragment.hasCentromere)						// If the fragment has centromere
        {
            centromereStrandID = fragment.oldStrandID;				// Store the associated original strand ID of that fragment
            ++centromereCount;							// Increment the centromere count for that data record
        }
    }

    if (centromereCount == 1)							// Only supports single-centromere records for determining home strand
    {										// Otherwise, home strand is the first listed fragment old strand ID
        return centromereStrandID;
    }

    // Zero or 2+ centromeres - ambiguous under today's single-slot
    // model. Fall back to the original strand-ordering convention.
    return record.fragments.empty() ? -1 : record.fragments[0].oldStrandID;
}



















// Mutated strands are only those with new strand IDs > numOriginalStrands.
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










// Finds record pairs representing a dicentric+acentric outcome, now
// generalized to records with any number of fragments (not just
// exactly 2), covering translocations exchanging segments from the
// middle of a chromosome, not just its ends.
//
// Within the dicentric candidate record, exactly two fragments must
// carry a centromere, from two DIFFERENT original strands, every
// other fragment in that record must be acentric and reference one
// of those same two strands. The first listed centromere-bearing
// fragment's strand becomes the dicentric record's own home strand;
// the second becomes the paired acentric record's home strand. This
// ordering is deliberate, whichever strand's fragment is listed first
// is the one that "owns" the dicentric drawing.
inline std::map<int, int> findDicentricAcentricHomeOverrides(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader)
{
    std::map<int, int> overrides;
    std::set<int> claimedAcentricRecords;				// Sorts Acentric IDs in order
    std::set<int> claimedExcisedRecords;				// Sorts excised fragment IDs and stores them

    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;
    const double tolerance = 0.001;					// In Mbp, amount fo discrepancy allowed between fragment ends and gap positions

    // Given pooled fragments for ONE strand (from a dicentric+acentric
    // pair), checks whether they (plus zero or more ADDITIONAL
    // single-fragment excised records for that strand) exactly span
    // its full declared length. This is what lets a deletion coexist
    // with the translocation that created the dicentric/acentric pair.
    // Only single-fragment excised pieces are matched (one deletion
    // per gap), a multi-fragment excised piece with its own internal
    // gaps isn't supported by this yet.
    auto tryTileWithDeletions = [&](std::vector<SDRfragment> pooled, int strandID,
                                     std::set<std::size_t> excludeIndices,
                                     std::vector<std::size_t>& usedExcisedIndices) -> bool
    {
        const std::size_t sizeIndex = static_cast<std::size_t>(strandID);
        if (sizeIndex >= masterHeader.intactChromosomeSizes.size())			// Ignore mutated chromosome entries, checking if intact chromosomes are spanned fully by the mutated segments within tolerance
        {
            return false;
        }

        const double fullSize = masterHeader.intactChromosomeSizes[sizeIndex];		// Size of chromosome in Mbp
        if (fullSize <= 0.0)
        {
            return false;
        }

	// Sort all 'pooled' fragments in a dicentric entry according to their start and end positions. First listed fragment with a centromere becomes
	// the home slot by which the dicentric chromosome will be drawn.
        std::sort(pooled.begin(), pooled.end(), [](const SDRfragment& a, const SDRfragment& b)
        {
            return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
        });


	// Find deletions that may affect clean tiling of fragments to the original chromosome's intact sizes
        struct Gap
        {
            double startPos;
            double endPos;
        };

        std::vector<Gap> gaps;

        const double firstStart = std::min(pooled.front().oldStartPosition, pooled.front().oldEndPosition);

        if (firstStart > 0.0 && !approxEqual(firstStart, 0.0, tolerance))
        {
            gaps.push_back({0.0, firstStart});
        }

        for (std::size_t k = 0; k + 1 < pooled.size(); ++k)
        {
            const double higherPos = std::max(pooled[k].oldStartPosition, pooled[k].oldEndPosition);
            const double lowerPos = std::min(pooled[k + 1].oldStartPosition, pooled[k + 1].oldEndPosition);

            if (approxEqual(higherPos, lowerPos, tolerance))
            {
                continue;
            }

            if (lowerPos <= higherPos)
            {
                return false;
            }

            gaps.push_back({higherPos, lowerPos});
        }

        const double lastEnd = std::max(pooled.back().oldStartPosition, pooled.back().oldEndPosition);

        if (fullSize > lastEnd && !approxEqual(lastEnd, fullSize, tolerance))
        {
            gaps.push_back({lastEnd, fullSize});
        }

        if (gaps.empty())				// Clean tiling already, no deletion involved.
        {
            return true;
        }

        std::vector<bool> gapMatched(gaps.size(), false);

        for (std::size_t r = 0; r < records.size(); ++r)
        {
            if (excludeIndices.count(r))
            {
                continue;
            }

            const SDRdataRecord& candidate = records[r];

            if (!candidate.linear || !isRearrangementCandidate(candidate.newStrandID, numOriginalStrands))
            {
                continue;
            }

            if (candidate.fragments.size() != 1)
            {
                continue;
            }

            const SDRfragment& fragment = candidate.fragments[0];

            if (fragment.oldStrandID != strandID || fragment.oldStartPosition > fragment.oldEndPosition)
            {
                continue;
            }

            for (std::size_t g = 0; g < gaps.size(); ++g)
            {
                if (gapMatched[g])
                {
                    continue;
                }

                if (approxEqual(fragment.oldStartPosition, gaps[g].startPos, tolerance) &&
                    approxEqual(fragment.oldEndPosition, gaps[g].endPos, tolerance))
                {
                    gapMatched[g] = true;
                    usedExcisedIndices.push_back(r);
                    break;
                }
            }
        }

        for (bool matched : gapMatched)
        {
            if (!matched)
            {
                return false;
            }
        }

        return true;
    };


    for (std::size_t i = 0; i < records.size(); ++i)
    {
        const SDRdataRecord& dicentricCandidate = records[i];

        if (!dicentricCandidate.linear || !isRearrangementCandidate(dicentricCandidate.newStrandID, numOriginalStrands))
        {
            continue;
        }

        std::vector<int> centromereStrandsInOrder;
        std::set<int> allStrandsInRecord;

        for (const SDRfragment& fragment : dicentricCandidate.fragments)
        {
            allStrandsInRecord.insert(fragment.oldStrandID);

            if (fragment.hasCentromere)
            {
                centromereStrandsInOrder.push_back(fragment.oldStrandID);
            }
        }

        if (centromereStrandsInOrder.size() != 2)
        {
            continue;
        }

        const int strandX = centromereStrandsInOrder[0];
        const int strandY = centromereStrandsInOrder[1];

        if (strandX == strandY || allStrandsInRecord.size() != 2)
        {
            continue;
        }

        for (std::size_t j = 0; j < records.size(); ++j)
        {
            if (i == j)
            {
                continue;
            }

            const SDRdataRecord& acentricCandidate = records[j];

            if (claimedAcentricRecords.count(acentricCandidate.newStrandID))
            {
                continue;
            }

            if (!acentricCandidate.linear || !isRearrangementCandidate(acentricCandidate.newStrandID, numOriginalStrands))
            {
                continue;
            }

            std::set<int> acentricStrands;
            bool hasCentromereFlag = false;

            for (const SDRfragment& fragment : acentricCandidate.fragments)
            {
                acentricStrands.insert(fragment.oldStrandID);

                if (fragment.hasCentromere)
                {
                    hasCentromereFlag = true;
                    break;
                }
            }

            if (hasCentromereFlag || acentricStrands.size() != 2 || !acentricStrands.count(strandX) || !acentricStrands.count(strandY))
            {
                continue;
            }

            std::vector<SDRfragment> pooledX;
            std::vector<SDRfragment> pooledY;

            for (const SDRfragment& fragment : dicentricCandidate.fragments)
            {
                (fragment.oldStrandID == strandX ? pooledX : pooledY).push_back(fragment);
            }

            for (const SDRfragment& fragment : acentricCandidate.fragments)
            {
                (fragment.oldStrandID == strandX ? pooledX : pooledY).push_back(fragment);
            }

            std::set<std::size_t> excludeIndices = {i, j};

            for (std::size_t r = 0; r < records.size(); ++r)
            {
                if (claimedExcisedRecords.count(records[r].newStrandID))
                {
                    excludeIndices.insert(r);
                }
            }

            std::vector<std::size_t> usedExcisedX;

            if (!tryTileWithDeletions(pooledX, strandX, excludeIndices, usedExcisedX))
            {
                continue;
            }

            std::set<std::size_t> excludeForY = excludeIndices;

            for (std::size_t idx : usedExcisedX)
            {
                excludeForY.insert(idx);
            }

            std::vector<std::size_t> usedExcisedY;

            if (!tryTileWithDeletions(pooledY, strandY, excludeForY, usedExcisedY))
            {
                continue;
            }

            overrides[dicentricCandidate.newStrandID] = strandX;
            overrides[acentricCandidate.newStrandID] = strandY;
            claimedAcentricRecords.insert(acentricCandidate.newStrandID);

            for (std::size_t idx : usedExcisedX)
            {
                claimedExcisedRecords.insert(records[idx].newStrandID);
            }

            for (std::size_t idx : usedExcisedY)
            {
                claimedExcisedRecords.insert(records[idx].newStrandID);
            }

            break;
        }
    }

    return overrides;
}











// Detects the 3-strand chromoplexy shape that findDicentricAcentricHomeOverrides()
// does NOT cover: strand A's centromeric piece fuses with strand B's centromeric
// piece (dicentric), but instead of strand B's OTHER piece pairing back up with
// strand A (which is what findDicentricAcentricHomeOverrides looks for), it
// instead fuses onto a THIRD strand C's centromeric piece. The result is
// monocentric overall (A keeps its own centromere via the dicentric record, C
// keeps its own via the second record), and strand B has been fully consumed,
// both of its pieces now live on strands A and C, with nothing left representing
// B on its own.
//
// Strand A's and strand C's own remaining acentric fragments (whatever did NOT
// fuse with B) are the two "loose ends" of this same rearrangement and are
// re-paired with each other here, to be drawn together in their own slot.
inline std::vector<SDRchromoplexyAcentricPlacement> findChromoplexyAcentricRecombinations(
    const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader)
{
    std::vector<SDRchromoplexyAcentricPlacement> placements;
    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;
    const double tolerance = 0.001;

    // Checks that fragment1 and fragment2 together exactly tile strandID's
    // whole declared length with no gap and no overlap - i.e. that strandID
    // really has been fully consumed by exactly these two pieces.
    auto fragmentsSpanFullLength = [&](const SDRfragment& fragment1, const SDRfragment& fragment2, int strandID) -> bool
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

        double lowerStart = std::min(fragment1.oldStartPosition, fragment1.oldEndPosition);
        double lowerEnd = std::max(fragment1.oldStartPosition, fragment1.oldEndPosition);
        double higherStart = std::min(fragment2.oldStartPosition, fragment2.oldEndPosition);
        double higherEnd = std::max(fragment2.oldStartPosition, fragment2.oldEndPosition);

        if (lowerStart > higherStart)
        {
            std::swap(lowerStart, higherStart);
            std::swap(lowerEnd, higherEnd);
        }

        return approxEqual(lowerStart, 0.0, tolerance)
            && approxEqual(lowerEnd, higherStart, tolerance)
            && approxEqual(higherEnd, fullSize, tolerance);
    };

    for (std::size_t i = 0; i < records.size(); ++i)
    {
        const SDRdataRecord& dicentric = records[i];

        if (!dicentric.linear || !isRearrangementCandidate(dicentric.newStrandID, numOriginalStrands) || dicentric.fragments.size() != 2)
        {
            continue;
        }

        std::vector<int> centromereStrands;

        for (const SDRfragment& fragment : dicentric.fragments)
        {
            if (fragment.hasCentromere)
            {
                centromereStrands.push_back(fragment.oldStrandID);
            }
        }

        if (centromereStrands.size() != 2 || centromereStrands[0] == centromereStrands[1])
        {
            continue; // Not a two-different-strand dicentric fusion.
        }

        const int strandA = centromereStrands[0];
        const int strandB = centromereStrands[1];

        const SDRfragment& fragmentA = (dicentric.fragments[0].oldStrandID == strandA) ? dicentric.fragments[0] : dicentric.fragments[1];
        const SDRfragment& fragmentB = (dicentric.fragments[0].oldStrandID == strandB) ? dicentric.fragments[0] : dicentric.fragments[1];

        // Look for a second record fusing strand B's OTHER (acentric) piece onto a third strand C's centromere.
        for (std::size_t j = 0; j < records.size(); ++j)
        {
            if (i == j)
            {
                continue;
            }

            const SDRdataRecord& monocentric = records[j];

            if (!monocentric.linear || !isRearrangementCandidate(monocentric.newStrandID, numOriginalStrands) || monocentric.fragments.size() != 2)
            {
                continue;
            }

            int strandC = -1;
            int centromereCountHere = 0;
            const SDRfragment* fragmentBOther = nullptr;
            const SDRfragment* fragmentC = nullptr;

            for (const SDRfragment& fragment : monocentric.fragments)
            {
                if (fragment.hasCentromere)
                {
                    ++centromereCountHere;

                    if (fragment.oldStrandID != strandB)
                    {
                        strandC = fragment.oldStrandID;
                        fragmentC = &fragment;
                    }
                }

                if (fragment.oldStrandID == strandB)
                {
                    fragmentBOther = &fragment;
                }
            }

            // Must be monocentric overall (strand C keeps exactly one centromere), must actually
            // reference strand B, that B fragment must be the acentric half, and C must be a genuinely
            // third strand, not A or B again.
            if (centromereCountHere != 1 || strandC == -1 || strandC == strandA || strandC == strandB
                || fragmentBOther == nullptr || fragmentBOther->hasCentromere || fragmentC == nullptr)
            {
                continue;
            }

            // Confirm strand B really is fully consumed: its two pieces (one from each record) tile its whole length.
            if (!fragmentsSpanFullLength(fragmentB, *fragmentBOther, strandB))
            {
                continue;
            }

            // Find strand A's and strand C's own leftover acentric fragments, whatever did NOT fuse
            // with strand B, and confirm each, together with its partner's centromeric fragment from
            // above, fully tiles its own strand's length. The SDR file may express these two leftover
            // pieces either as ONE record that already combines both (fragments.size() == 2, one
            // fragment per strand) or as two separate single-fragment records, one per strand.
            std::vector<int> leftoverNewStrandIDs;

            for (const SDRdataRecord& candidate : records)
            {
                if (candidate.newStrandID == dicentric.newStrandID || candidate.newStrandID == monocentric.newStrandID)
                {
                    continue;
                }

                if (!candidate.linear || !isRearrangementCandidate(candidate.newStrandID, numOriginalStrands))
                {
                    continue;
                }

                if (candidate.fragments.size() == 2)
                {
                    // Already-combined shape: one fragment must be strand A's leftover, the other
                    // strand C's leftover, neither carrying a centromere.
                    const SDRfragment* candidateA = nullptr;
                    const SDRfragment* candidateC = nullptr;

                    for (const SDRfragment& fragment : candidate.fragments)
                    {
                        if (fragment.hasCentromere)
                        {
                            candidateA = nullptr;
                            candidateC = nullptr;
                            break;
                        }

                        if (fragment.oldStrandID == strandA)
                        {
                            candidateA = &fragment;
                        }
                        else if (fragment.oldStrandID == strandC)
                        {
                            candidateC = &fragment;
                        }
                    }

                    if (candidateA != nullptr && candidateC != nullptr
                        && fragmentsSpanFullLength(fragmentA, *candidateA, strandA)
                        && fragmentsSpanFullLength(*fragmentC, *candidateC, strandC))
                    {
                        leftoverNewStrandIDs = { candidate.newStrandID };
                        break;
                    }
                }
            }

            if (leftoverNewStrandIDs.empty())
            {
                // Fall back to the two-separate-records shape.
                const SDRdataRecord* leftoverA = nullptr;
                const SDRdataRecord* leftoverC = nullptr;

                for (const SDRdataRecord& candidate : records)
                {
                    if (candidate.newStrandID == dicentric.newStrandID || candidate.newStrandID == monocentric.newStrandID)
                    {
                        continue;
                    }

                    if (!candidate.linear || !isRearrangementCandidate(candidate.newStrandID, numOriginalStrands) || candidate.fragments.size() != 1)
                    {
                        continue;
                    }

                    const SDRfragment& fragment = candidate.fragments[0];

                    if (fragment.hasCentromere)
                    {
                        continue;
                    }

                    if (leftoverA == nullptr && fragment.oldStrandID == strandA && fragmentsSpanFullLength(fragmentA, fragment, strandA))
                    {
                        leftoverA = &candidate;
                    }
                    else if (leftoverC == nullptr && fragment.oldStrandID == strandC && fragmentsSpanFullLength(*fragmentC, fragment, strandC))
                    {
                        leftoverC = &candidate;
                    }
                }

                if (leftoverA != nullptr && leftoverC != nullptr)
                {
                    leftoverNewStrandIDs = { leftoverA->newStrandID, leftoverC->newStrandID };
                }
            }

            if (leftoverNewStrandIDs.empty())
            {
                continue;
            }

            SDRchromoplexyAcentricPlacement placement{};
            placement.consumedStrandID = strandB;
            placement.leftoverNewStrandIDs = leftoverNewStrandIDs;
            placements.push_back(placement);

            break; // Found this dicentric record's match, move on to the next candidate dicentric record.
        }
    }

    return placements;
}


















// --------------------------------------------------------------------------- //
// Detects LONG DELETION events within a single cell.
// --------------------------------------------------------------------------- //

// Modified detectDeletions function to detect multiple deletions within a single strand. A strand with N deletions will contain N + 1 fragments in SDR data field 3, and will have
// N + 1 data entries, 1 entry being the original strand with all the gaps, and N entries representing each deletion causing the gaps.
inline std::vector<SDRdeletionEvent> detectDeletions(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader)
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

	bool multipleRemainingRecords = false;					// Variable to check if the amount of fragments remaining in a strand that
										// had an excised fragment is greater than 1.

        for (auto& [newStrandID, fragments] : groupNewStrand)			// Loop through each pair in groupNewStrand vector
        {
	    if (fragments.size() != newStrandIDtotalFragments[newStrandID]) 	// Reject if this record has fragments from OTHER old strands too - a pure remaining/excised piece must be entirely from this old strand.
            {
                continue;
            }

            if (fragments.size() >= 2)						// If the entry has two flanking fragments, then there are multiple
            {									// remaining records. A deletion was detected.
		if (flankingFragments != nullptr)
		{
		    multipleRemainingRecords = true;
		    break;
		}

		remainingStrandID = newStrandID;                                // Track the newStrandID of the flanking fragments in the remaining strand
                flankingFragments = &fragments;
            }
            else if (fragments.size() == 1)					// Fragment size of 1 corresponds to the excised/deleted DNA fragment
            {
		excisedCandidates.push_back({newStrandID, &fragments[0]});	// Store excised candidate fragment new Strand ID and corresponding old Strand ID.

            }
        }

	if (multipleRemainingRecords || flankingFragments == nullptr || excisedCandidates.empty())
        {
            continue;
        }


        // Sort the flanking fragments by start position so we can find the gap between them regardless of file order.
        std::sort(flankingFragments->begin(), flankingFragments->end(),[](const SDRfragment& a, const SDRfragment& b)
	    {
                return a.oldStartPosition < b.oldStartPosition;
            });


	// -------------------------------------------------------- //
	// Check if fragments are missing at the ends or the center
	// -------------------------------------------------------- //

	// The chromosome's true declared length is needed to detect a
        // TERMINAL deletion, one touching position 0 or the chromosome's
        // very end. A terminal deletion leaves no extra "remaining"
        // fragment on that side (there's nothing before position 0, or
        // after the chromosome's end, to remain). So the old rigid
        // "N+1 remaining fragments for N deletions" assumption only
        // held when every deletion was strictly internal.
        const std::size_t sizeIndex = static_cast<std::size_t>(oldStrandID);
        const double chromosomeFullSizeMbp = (sizeIndex < masterHeader.intactChromosomeSizes.size()) ? masterHeader.intactChromosomeSizes[sizeIndex] : 0.0;

	struct Gap								// Struct to store the size of the gap from the excised piece
	{									// and match it with the flanking fragment ends within tolerance
	    double startPos;
	    double endPos;
	};

	std::vector<Gap> gaps;
	bool validGaps = true;

	// Leading gap from the chromosome's true start (0) to the
        // first flanking fragment's start.
        if (chromosomeFullSizeMbp > 0.0 && !approxEqual((*flankingFragments)[0].oldStartPosition, 0.0, delTolerance))
        {
            const double gapStart = 0.0;
            const double gapEnd = (*flankingFragments)[0].oldStartPosition;

            if (gapEnd > gapStart)
            {
                gaps.push_back({gapStart, gapEnd});
            }
        }


	// Internal gaps between consecutive flanking fragments.
	for (std::size_t i = 0; i + 1 < flankingFragments->size(); i++)
        {
	    // Generalized for multiple flanking fragments per data record
            const double gapStart = (*flankingFragments)[i].oldEndPosition;	// Gap starts at the end of the first flanking fragment
            const double gapEnd = (*flankingFragments)[i + 1].oldStartPosition; // Gap ends at the start of the second flanking fragment

            if (gapEnd <= gapStart)                                             // Check if a gap was left from the deletion event
            {
                validGaps = false;
                break;                                                          // Overlapping/non-gap boundary - not a valid deletion shape for this strand.
            }

            gaps.push_back({gapStart, gapEnd});
        }

	if (!validGaps)
	{
	    continue;
	}



	// Trailing gap from the last flanking fragment's end to the
        // chromosome's true full length.
        if (chromosomeFullSizeMbp > 0.0)
        {
            const double lastFragmentEnd = flankingFragments->back().oldEndPosition;

            if (chromosomeFullSizeMbp > lastFragmentEnd && !approxEqual(lastFragmentEnd, chromosomeFullSizeMbp, delTolerance))
            {
                gaps.push_back({lastFragmentEnd, chromosomeFullSizeMbp});
            }
        }


	if (gaps.size() != excisedCandidates.size())
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

    const double balTraTolerance = 0.00001;						// Tolerance of base position mismatch of approximately 10 bp


    // Sub function to check if the strand being investigated is valid (linear for BalTrans, and a mutation record). Need this since we are
    // checking across multiple data records for two distinct old strand IDs sharing DNA information in a balanced manner (no information loss).
    // Groups a given record's fragments by each old Strand ID that fragment came from
    auto getRecordStrandFragments = [&](const SDRdataRecord& record, std::map<int, std::vector<SDRfragment>>& groupStrand) -> bool
    {
        if (!record.linear)								// Data record must be linear to be considered
        {
            return false;
        }

        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))		// Data record must have a new strand ID > numOriginalStrands
        {										// Otherwise not a mutation
            return false;
        }

        for (const SDRfragment& fragment : record.fragments)				// If initial checks pass, store the fragment associated with
        {										// A given oldStrandID
            groupStrand[fragment.oldStrandID].push_back(fragment);
        }

        return !groupStrand.empty();							// Return true if the groupStrand is non-empty
    };



    // Checks that the pooled fragments for one strand (groupStrand), normalized to [min,max] and sorted, span the full size of 
    // the chromosome [0, fullSize] exactly with no gaps or overlaps. On success, breakpointsOut holds every internal
    // boundary between consecutive fragments.
    auto checkChromosomeSpan = [&](std::vector<SDRfragment> pooledFragments, int strandID, std::vector<double>& allBreakpoints) -> bool
    {
        const std::size_t sizeIndex = static_cast<std::size_t>(strandID);

        if (sizeIndex >= masterHeader.intactChromosomeSizes.size())			// Check sizes of intact chromosomes to ensure translocations are
        {										// balanced, will skip all mutated data entries
            return false;
        }

        const double fullSize = masterHeader.intactChromosomeSizes[sizeIndex];		// Store the fullSize of a given chromosome

        if (fullSize <= 0.0)
        {
            return false;
	}


	struct NormalizedFragment							// Struct to order fragments according to start position
        {
            double lowerPos;								// Lower of the two fragment positions (generalized in case fragment inverted)
            double higherPos;								// Higher of the two fragment positions
        };

        std::vector<NormalizedFragment> normalizedFragments;				// Vector to hold the rearranged (normalized) fragment order after
        normalizedFragments.reserve(pooledFragments.size());

	for (const SDRfragment& fragment : pooledFragments)				// Loop through all pooledFragments which group fragments according to old strand ID
        {
            normalizedFragments.push_back({						// Store the fragments according to lower fragment position first,
                std::min(fragment.oldStartPosition, fragment.oldEndPosition),		// higher fragment position second.
                std::max(fragment.oldStartPosition, fragment.oldEndPosition)
            });
        }

	// Sort normalizedFragments vector fragments in ascending order of position.
        std::sort(normalizedFragments.begin(), normalizedFragments.end(), [](const NormalizedFragment& a, const NormalizedFragment& b)
        {
            return a.lowerPos < b.lowerPos;
        });

	// Need the first sorted fragment's start position to start at 0
        if (!approxEqual(normalizedFragments.front().lowerPos, 0.0, balTraTolerance))
        {
            return false;
        }

	// Need the last sorted fragment's end position to end at the fullSize of the chromosome
	if (!approxEqual(normalizedFragments.back().higherPos, fullSize, balTraTolerance))
        {
            return false;
        }


        allBreakpoints.clear();

	// Loop through all normalizedFragments to check if the fragments are contiguous
        for (std::size_t k = 0; k + 1 < normalizedFragments.size(); ++k)
        {
	    // Require neighbouring fragments to be contiguous within balTraTolerance = 10 bp, otherwise not a balanced translocation
            if (!approxEqual(normalizedFragments[k].higherPos, normalizedFragments[k + 1].lowerPos, balTraTolerance))
            {
                return false;
            }
	    // Store breakpoint locations using the higher positions of the fragment ends
            allBreakpoints.push_back(normalizedFragments[k].higherPos);
        }

        return true;
    };



    for (std::size_t i = 0; i < records.size(); ++i)				// Loop through all SDR data records, looking for a first candidate
    {
        std::map<int, std::vector<SDRfragment>> firstStrandGroup;		// This variable stores the first group of two original strands contributing to the 
										// balanced translocation.
        if (!getRecordStrandFragments(records[i], firstStrandGroup))		// If the first record does not fit the required balTrans shape, skip record
        {
            continue;
        }

        if (firstStrandGroup.size() != 2)					// Simple balanced translocation mutations only involve two chromosomes
        {
            continue;
        }

	for (std::size_t j = i + 1; j < records.size(); ++j)			// Loop through all following records after the i-th data record,
        {									// looking for a second candidate for the balTrans
            std::map<int, std::vector<SDRfragment>> secondStrandGroup;		// Store the second group of two strands that contribute to balanced translocations.

            if (!getRecordStrandFragments(records[j], secondStrandGroup))	// If secondStrandGroup does not fit balTrans shape, skip record
            {
                continue;
            }

            if (secondStrandGroup.size() != 2)					// Simple balTrans mutations only involve two chromosomes
            {
                continue;
            }


	    // Store fragment original strand IDs for the first candidate
            std::vector<int> firstStrands;

	    for (auto& [strandID, fragments] : firstStrandGroup)
            {
                firstStrands.push_back(strandID);
            }

	    // Store the fragment original strand IDs for the second candidate
            std::vector<int> secondStrands;

            for (auto& [strandID, fragments] : secondStrandGroup)
            {
                secondStrands.push_back(strandID);
            }

	    // If the two original strands referenced by each candidate record are not equivalent, then this is not a balanced translocation
            if (firstStrands != secondStrands)
            {
                continue;
            }

	    // Once the strand IDs are confirmed to be identical, set strand A and strand B IDs
	    const int strandAid = firstStrands[0];			// equivalent to secondStrand[0]
            const int strandBid = firstStrands[1];			// equivalent to secondStrand[1]


	    // Instead of assuming each record has exactly one fragment per strand (which only holds for a clean single-breakpoint swap), 
	    // this gathers every fragment referencing strand A across both records into one list, and does the same for strand B.
	    // This is what also allows for central segment balanced translocations, not only chromosome end swapping.
            std::vector<SDRfragment> pooledFragmentsA = firstStrandGroup[strandAid];
            pooledFragmentsA.insert(pooledFragmentsA.end(), secondStrandGroup[strandAid].begin(), secondStrandGroup[strandAid].end());

            std::vector<SDRfragment> pooledFragmentsB = firstStrandGroup[strandBid];
            pooledFragmentsB.insert(pooledFragmentsB.end(), secondStrandGroup[strandBid].begin(), secondStrandGroup[strandBid].end());

	    // Store the breakpoints on both chromosome A and chromosome B
            std::vector<double> breakpointsA;
            std::vector<double> breakpointsB;

	    // If the fragments of a given old strand ID do not span the entire length of the original intact chromosome they came from,
	    // this is not a balanced translocation
            if (!checkChromosomeSpan(pooledFragmentsA, strandAid, breakpointsA) || !checkChromosomeSpan(pooledFragmentsB, strandBid, breakpointsB))
            {
                continue;
            }


	    // All balanced translocation mutation checks passed, store the strand and fragment information, as well as the break locations
	    SDRtranslocationEvent event{};
            event.oldStrandAid = strandAid;
            event.oldStrandBid = strandBid;
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

inline std::vector<SDRecDNAevent> detectECDNA(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader)
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
        std::vector<std::pair<int, std::vector<SDRfragment>*>> excisedCandidates; // {newStrandID, fragments}


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
                    multipleRemainingRecords = true; 				// More than one candidate "remaining" record - invalid shape.
                    break;
                }

                remainingStrandID = newStrandID;
                flankingFragments = &fragments;
            }
            else if (!record->linear)						// Necessary format is one entry that is non-linear with one singular fragment referencing the same oldStrandID as the record with the two flanking fragments
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


	// ----------------------------------------------------------- //
	// Determine if the deleted segment is from the center or ends of the chromosome
	// ----------------------------------------------------------- //

	// A TERMINAL ecDNA excision (touching position 0 or the
        // chromosome's true end) leaves no extra flanking fragment on
        // that side, same as a terminal long deletion, checking
        // against the chromosome's true declared length generalizes
        // this beyond just gaps between listed fragments.

	// Checks if deleted fragment is at the start of the chromosome
        const std::size_t sizeIndex = static_cast<std::size_t>(oldStrandID);
        const double chromosomeFullSizeMbp = (sizeIndex < masterHeader.intactChromosomeSizes.size()) ? masterHeader.intactChromosomeSizes[sizeIndex] : 0.0;

        if (chromosomeFullSizeMbp > 0.0 && !approxEqual(flankingFragments->front().oldStartPosition, 0.0, ecDNAtolerance)
            && flankingFragments->front().oldStartPosition > 0.0)
        {
            gaps.push_back({0.0, flankingFragments->front().oldStartPosition});
        }

	// Checks if deleted segment is in the central section of the chromosome.
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


	// Checks if deleted segment is at the end of the chromosome
	if (chromosomeFullSizeMbp > 0.0)
        {
            const double lastFragmentEnd = flankingFragments->back().oldEndPosition;

            if (chromosomeFullSizeMbp > lastFragmentEnd && !approxEqual(lastFragmentEnd, chromosomeFullSizeMbp, ecDNAtolerance))
            {
                gaps.push_back({lastFragmentEnd, chromosomeFullSizeMbp});
            }
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

// Two candidate records, each with exactly 2 fragments from 2 different old strands
// (linear, any orientation, a reversed fragment represents an inverted exchange, e.g.
// one forming a dicentric/acentric pair, and is not itself a reason to reject). For each
// of the two strands involved, the two records' fragments for that strand are compared
// using NORMALIZED [min,max] bounds (so a reversed fragment's boundaries are still
// checked correctly, unlike a raw oldStart/oldEnd comparison): if they touch exactly,
// that strand is "clean" (no deletion); if BOTH strands are clean, this is a plain
// balanced translocation, handled elsewhere. If exactly one strand has a genuine gap
// instead, a third record (single fragment, from the gapped strand, non-reversed,
// exactly filling the gap) confirms a deletion-translocation.

inline std::vector<SDRdeletionTranslocationEvent> detectDeletionTranslocations(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader)
{
    std::vector<SDRdeletionTranslocationEvent> delTras;					// Store deletion-translocation events

    const std::vector<SDRdataRecord>& records = subHeader.dataRecords;
    const double tolerance = 0.001;							// Deletion base mismatch tolerance ~ 1000 bp = 0.001 Mbp

    // Single-use function to determine which fragments in a record should be considered for deletion-translocation.
    auto getRecordStrandFragments = [&](const SDRdataRecord& record, std::map<int, std::vector<SDRfragment>>& byStrand) -> bool
    {
        if (!record.linear || !isRearrangementCandidate(record.newStrandID, numOriginalStrands))
        {
            return false;
        }

        for (const SDRfragment& fragment : record.fragments)
        {
            byStrand[fragment.oldStrandID].push_back(fragment);
        }

        return !byStrand.empty();
    };

    // Checks whether the pooled fragments for one strand (normalized,
    // sorted) tile [0, fullSize] exactly "clean". If there's exactly
    // ONE gap among them (and the pooled fragments still start at 0
    // and end at fullSize), returns false with hasGap set, the
    // deletion candidate for that strand. Any other shape (doesn't
    // start at 0, more than one gap, an overlap) is rejected outright.
    auto checkStrandTiling = [&](std::vector<SDRfragment> pooled, int strandID,
                              std::vector<double>& breakpointsOut,
                              bool& hasGap, double& gapStart, double& gapEnd) -> bool
    {
    	breakpointsOut.clear();
    	hasGap = false;

    	const std::size_t sizeIndex = static_cast<std::size_t>(strandID);

    	if (sizeIndex >= masterHeader.intactChromosomeSizes.size())			// Skip intact chromosome entries
    	{
            return false;
    	}

    	const double fullSize = masterHeader.intactChromosomeSizes[sizeIndex];		// The size of each chromosome in Mbp

    	if (fullSize <= 0.0)
    	{
            return false;
    	}

	// Sort the fragments in a given entry in ascending order of fragment start and end position (normalized for inverted strands)
    	std::sort(pooled.begin(), pooled.end(), [](const SDRfragment& a, const SDRfragment& b)
    	{
            return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
    	});

    	int gapCount = 0;

    	// Leading gap: before the first fragment (a deletion touching the chromosome's very start).
    	const double firstGapLower = std::min(pooled.front().oldStartPosition, pooled.front().oldEndPosition);

    	if (firstGapLower > 0.0 && !approxEqual(firstGapLower, 0.0, tolerance))
    	{
            ++gapCount;
            hasGap = true;
            gapStart = 0.0;
            gapEnd = firstGapLower;
    	}

    	// Internal gaps: between consecutive fragments.
    	for (std::size_t k = 0; k + 1 < pooled.size(); ++k)
    	{
            const double higherFragmentPos = std::max(pooled[k].oldStartPosition, pooled[k].oldEndPosition);
            const double lowerFragmentPos = std::min(pooled[k + 1].oldStartPosition, pooled[k + 1].oldEndPosition);

            if (approxEqual(higherFragmentPos, lowerFragmentPos, tolerance))
            {
            	breakpointsOut.push_back(higherFragmentPos);
            	continue;
            }

            if (lowerFragmentPos <= higherFragmentPos)		// Overlap = invalid shape.
            {
            	return false;
            }

            ++gapCount;

            if (gapCount > 1)					// Only one deletion per strand is supported here.
            {
            	return false;
            }

            hasGap = true;
            gapStart = higherFragmentPos;
            gapEnd = lowerFragmentPos;
    	}

    	// Trailing gap: after the last fragment (a deletion touching the
    	// chromosome's very end).
    	const double lastGapHigher = std::max(pooled.back().oldStartPosition, pooled.back().oldEndPosition);

    	if (fullSize > lastGapHigher && !approxEqual(lastGapHigher, fullSize, tolerance))
    	{
            ++gapCount;

            if (gapCount > 1)
            {
            	return false;
            }

            hasGap = true;
            gapStart = lastGapHigher;
            gapEnd = fullSize;
    	}

	// true = clean tiling (no gaps anywhere); false (with hasGap set) = the single-deletion candidate shape, wherever the gap fell.
    	return !hasGap;
    };


    for (std::size_t i = 0; i < records.size(); ++i)					// Loop through all SDR data records
    {
        std::map<int, std::vector<SDRfragment>> firstStrandGroup;			// Look at the first mutated record ID and its fragments

	// Check if SDR record is a mutation and if it is composed of fragments from exactly two different original strands
        if (!getRecordStrandFragments(records[i], firstStrandGroup) || firstStrandGroup.size() != 2)
        {
            continue;
        }


        for (std::size_t j = i + 1; j < records.size(); ++j)				// Loop through all the following SDR data records after the i-th record
        {
            std::map<int, std::vector<SDRfragment>> secondStrandGroup;			// Look at the second mutated record ID and its fragments

	    // Check if SDR record is a mutation and if it is composed of fragments from exactly two different original strands
            if (!getRecordStrandFragments(records[j], secondStrandGroup) || secondStrandGroup.size() != 2)
            {
                continue;
            }

            std::vector<int> firstStrands;					// Vector storing the IDs of the fragments making up the first data record
            for (auto& [strandID, fragments] : firstStrandGroup) 
	    { 
		firstStrands.push_back(strandID); 
	    }

            std::vector<int> secondStrands;					// Vector storing the IDs of the fragments making up the second data record
            for (auto& [strandID, fragments] : secondStrandGroup) 
	    { 
		secondStrands.push_back(strandID); 
	    }

	    // If the old strand IDs in both records are not exaclty the same, this shape is invalid, skip this record
            if (firstStrands != secondStrands)
            {
                continue;
            }


            const int strandAid = firstStrands[0];	// Equal to secondStrands[0]
            const int strandBid = firstStrands[1];	// Equal to secondStrands[1]

	    // Pooled A stores all fragments with old strand ID 'strandAid' from both data records involved
            std::vector<SDRfragment> pooledA = firstStrandGroup[strandAid];
            pooledA.insert(pooledA.end(), secondStrandGroup[strandAid].begin(), secondStrandGroup[strandAid].end());
	    // Pooled B stores all fragments with old strand ID 'strandBid' from both data records involved
            std::vector<SDRfragment> pooledB = firstStrandGroup[strandBid];
            pooledB.insert(pooledB.end(), secondStrandGroup[strandBid].begin(), secondStrandGroup[strandBid].end());

            std::vector<double> breakpointsA;
            bool hasGapA = false;
            double gapStartA = 0.0, gapEndA = 0.0;
            const bool cleanA = checkStrandTiling(pooledA, strandAid, breakpointsA, hasGapA, gapStartA, gapEndA);

            std::vector<double> breakpointsB;
            bool hasGapB = false;
            double gapStartB = 0.0, gapEndB = 0.0;
            const bool cleanB = checkStrandTiling(pooledB, strandBid, breakpointsB, hasGapB, gapStartB, gapEndB);

            if (cleanA && cleanB)			// Plain balanced translocation, handled elsewhere.
            {
                continue;
            }

            double gapStart = 0.0;
	    double gapEnd = 0.0;
            int deletedOldStrandID = -1;
            std::vector<double> cleanBreakPositions;

	    // Check to make sure a deletion only exists on one of the two strands (strand A clean, strand B has deletion)
            if (cleanA && !cleanB && hasGapB)
            {
                gapStart = gapStartB;
                gapEnd = gapEndB;
                deletedOldStrandID = strandBid;
                cleanBreakPositions = breakpointsA;
            }
            else if (cleanB && !cleanA && hasGapA)	// (strand B clean, strand A has deletion)
            {
                gapStart = gapStartA;
                gapEnd = gapEndA;
                deletedOldStrandID = strandAid;
                cleanBreakPositions = breakpointsB;
            }
            else					// Cannot have a deletion in both strand A and strand B with a translocation, wrong shape
            {
                continue;
            }


	    // Checking geometry of mutation entries.
            for (const SDRdataRecord& deletionRecord : records)
            {
		// Do not count lone deletion records in the deletion-translocation checks
                if (deletionRecord.newStrandID == records[i].newStrandID || deletionRecord.newStrandID == records[j].newStrandID)
                {
                    continue;
                }
		// Deleted segments must be linear and be a mutated entry
                if (!deletionRecord.linear || !isRearrangementCandidate(deletionRecord.newStrandID, numOriginalStrands))
                {
                    continue;
                }
		// Resulting deletion can only be one fragment
                if (deletionRecord.fragments.size() != 1)
                {
                    continue;
                }

                const SDRfragment& excisedFragment = deletionRecord.fragments[0];
		// The deleted segment must come from the same original chromosome as the deletion-translocation
                if (excisedFragment.oldStrandID != deletedOldStrandID)
                {
                    continue;
                }

                const double excisedLowerPos = std::min(excisedFragment.oldStartPosition, excisedFragment.oldEndPosition);
                const double excisedHigherPos = std::max(excisedFragment.oldStartPosition, excisedFragment.oldEndPosition);

		// The excised fragment ends should match within tolerance the gap region
                if (!approxEqual(excisedLowerPos, gapStart, tolerance) || !approxEqual(excisedHigherPos, gapEnd, tolerance))
                {
                    continue;
                }

		// All check passed, store deletion-translocation event
                SDRdeletionTranslocationEvent event{};
                event.oldStrandAid = strandAid;
                event.oldStrandBid = strandBid;
                event.deletedOldStrandID = deletedOldStrandID;
                event.cleanBreakPositions = cleanBreakPositions;
                event.deletionStart = gapStart;
                event.deletionEnd = gapEnd;
                event.newStrandID1 = records[i].newStrandID;
                event.newStrandID2 = records[j].newStrandID;
                event.excisedStrandID = deletionRecord.newStrandID;

                delTras.push_back(event);
                break;
            }
        }
    }

    return delTras;
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
// PAIR of distinct old strand IDs referenced by its fragments gets an edge, regardless
// of whether that record matches any specific named mutation shape (translocation,
// deletion-translocation, etc.). Any connected component spanning 3 or more distinct
// strands is reported as one chromoplexy event.
inline std::vector<SDRchromoplexyEvent> detectChromoplexy(const SDRsubHeader& subHeader, int numOriginalStrands, const SDRmasterHeader& masterHeader)
{
    std::vector<SDRchromoplexyEvent> chromoplexyEvents;

    const std::vector<SDRdeletionEvent> deletions = detectDeletions(subHeader, numOriginalStrands, masterHeader);


    // Initial State: Every chromosome is isolated in its own group (parent[x] = x).
    // Mutation 1 (Chromosome 1 and 2 swap material): You call unite(1, 2). Chromosome 1 and 2 are now linked. parent[1] becomes 2.
    // Mutation 2 (Chromosome 2 and 3 swap material): You call unite(2, 3). Now, Chromosome 1, 2, and 3 are all part of the same group led by 3.
    // If you later ask findRoot(1), the code will climb up (1 -> 2 -> 3), optimize the path, and return 3.

    // Groups chromosomes based on their shared mutations
    std::map<int, int> parent;

    auto findRoot = [&](int x) -> int		// x are the individual chromosome IDs
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
    std::vector<std::pair<int, std::set<int>>> recordStrandSets; 			// {newStrandID, {distinct old strand IDs}}

    for (const SDRdataRecord& record : subHeader.dataRecords)				// Loop through all SDR data records per cell
    {
        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))		// If the record is not a mutation, but an intact record, skip
        {
            continue;
        }

	std::set<int> strandsInRecord;							// Track which chromosomes make up the fragments in the mutation

        for (const SDRfragment& fragment : record.fragments)
        {
            strandsInRecord.insert(fragment.oldStrandID);
        }

        if (strandsInRecord.size() < 2)							// Purely intra-strand record, contributes no edges (only one chromosome is involved).
        {
            continue;
        }

        recordStrandSets.push_back({record.newStrandID, strandsInRecord});		// Couple the new strand ID with the old strand IDs of its fragments

        std::vector<int> strandVec(strandsInRecord.begin(), strandsInRecord.end());	// Store the involved new strand IDs and their respective fragments' old strand IDs

	// Loop through all original chromosome IDs, in the current record and following records
	for (std::size_t a = 0; a < strandVec.size(); ++a)
        {
            for (std::size_t b = a + 1; b < strandVec.size(); ++b)
            {
                unite(strandVec[a], strandVec[b]);					// Unite two old strand IDs from two records
            }
        }
    }


    // as you iterate all strands, every strand sharing the same root ends up collected together in one set — and using a set rather than a vector means 
    // duplicates (from path-compressed lookups) just don't matter.
    std::map<int, std::set<int>> componentsByRoot;					// maps one arbitrary representative strand the full set of every strand connected to it.

    // for every strand that's ever been registered as a node: call findRoot(strandID) to get its component's representative 
    // ("root") strand ID, and record that strandID belongs under that root.
    for (auto& [strandID, strandParent] : parent)
    {
        componentsByRoot[findRoot(strandID)].insert(strandID);
    }

    // Loops through every discovered component and checks if it has >= 3 strands
    for (auto& [root, strandSet] : componentsByRoot)
    {
        if (strandSet.size() < 3)							// Chromoplexy involved 3 or more chromosomes
        {
            continue;
        }

	// Build the list of records by newStrandID since the component qualifies with 3+ strands
        std::vector<int> contributingIDs;						// Will hold all records touching at least one of the componentsByRoot's 3+ strands.

	// The nested loop checks for each record if any old strand this record touches belong to the current 3+-strand component.
        for (auto& [newStrandID, strands] : recordStrandSets)
        {
	    // Check each strand reference in the given record
            for (int s : strands)
            {
                if (strandSet.count(s))							// If that strand is in strandSet, add newStrandID to the contributingIDs
                {
                    contributingIDs.push_back(newStrandID);				// Store the record ID and break, one match is proof the record belongs to the chromoplexy event, stop searching further.
                    break;
                }
            }
        }


	// All checks passed for chromoplexy, store the mutation event information
        SDRchromoplexyEvent event{};
        event.involvedStrandIDs = std::vector<int>(strandSet.begin(), strandSet.end());
        event.contributingNewStrandIDs = contributingIDs;


	// Extra non-essential check to assess whether one of the involved strands in chromoplexy also has a separate associated deletion, and stores that deletions information
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
// rearranged record touching the cluster's strand(s), a direct measure of how many
// pieces the strand(s) were actually broken into, regardless of shape. Only counts the 
// remaining fragments in the original strand, not the excised pieces, as a
// clearer metric for the number of rearrangements that occurred in a strand.
inline std::vector<SDRchromothripsisEvent> detectChromothripsis(const SDRsubHeader& subHeader, int numOriginalStrands, int minFragmentCount = 10)
{
    std::vector<SDRchromothripsisEvent> chromothripsisEvents;			// Store chromothripsis events

    // Groups chromosomes based on their shared mutations
    std::map<int, int> parent;

    auto findRoot = [&](int x) -> int		// x are the individual chromosome IDs
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


    // Map from old strand ID -> set of newStrandIDs (records) whose
    // fragments touch it, whether intra-strand-only or inter-strand.
    std::map<int, std::set<int>> recordsByStrand;

    // Map from newStrandID -> the record's own total fragment count, so
    // each record's fragments are only ever tallied once even if it
    // touches multiple strands in a cluster.
    std::map<int, std::size_t> fragmentCountByNewStrandID;

    for (const SDRdataRecord& record : subHeader.dataRecords)
    {
	// If record is not a mutated record, skip
        if (!isRearrangementCandidate(record.newStrandID, numOriginalStrands))
        {
            continue;
        }

	// Old strand IDs of all fragments in a given record
        std::set<int> strandsInRecord;

        for (const SDRfragment& fragment : record.fragments)
        {
            strandsInRecord.insert(fragment.oldStrandID);
        }

	// Count how many fragments there are across all involved new Strand IDs in a chromothripsis event (either 1 or 2 new Strand ID, i.e. data records)
        fragmentCountByNewStrandID[record.newStrandID] = record.fragments.size();

        for (int strandID : strandsInRecord)
        {
            findRoot(strandID); 				// Register the node even if this record touches only one strand.
            recordsByStrand[strandID].insert(record.newStrandID);
        }

        if (strandsInRecord.size() < 2)				// Intra-strand record, no edge to add.
        {
            continue;
        }

	// Store all the strand IDs involved across all chromothripsis-involved records
        std::vector<int> strandVec(strandsInRecord.begin(), strandsInRecord.end());

	// Unite a and b data records that could involve two old strand IDs in a chromothripsis event
        for (std::size_t a = 0; a < strandVec.size(); a++)
        {
            for (std::size_t b = a + 1; b < strandVec.size(); b++)
            {
                unite(strandVec[a], strandVec[b]);
            }
        }
    }

    // Map the new Strand ID to the old Strand IDs of the fragments that make it up
    std::map<int, std::set<int>> componentsByRoot;

    for (auto& [strandID, strandParent] : parent)
    {
        componentsByRoot[findRoot(strandID)].insert(strandID);
    }


    for (auto& [root, strandSet] : componentsByRoot)
    {
        if (strandSet.size() > 2)
        {
            continue; // 3+ strands, chromoplexy's territory, not chromothripsis.
        }

        // Any record with 2+ fragments is a "major" surviving piece and
        // counts in full, whether it's a single-strand remaining piece
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
                if (fragmentCountByNewStrandID[newStrandID] >= 2)		// Ignores single excised fragments, only strands made up of multiple fragments count
                {								// towards chromothripsis
                    majorNewStrandIDs.insert(newStrandID);
                }
            }
        }

        int totalFragmentCount = 0;						// Count how many contributing fragments there are to be considered as chromothripsis

        for (int newStrandID : majorNewStrandIDs)				// Count all non-single-excised fragment pieces across the involved parent strands
        {
            totalFragmentCount += static_cast<int>(fragmentCountByNewStrandID[newStrandID]);
        }

	// Min fragment count = 10 by default, function argument parameter can be changed.
        if (totalFragmentCount < minFragmentCount)
        {
            continue;
        }

	// Checks passed, store the chromothripsis event
        SDRchromothripsisEvent event{};
        event.involvedStrandIDs = std::vector<int>(strandSet.begin(), strandSet.end());
        event.contributingNewStrandIDs = std::vector<int>(majorNewStrandIDs.begin(), majorNewStrandIDs.end());
        event.totalFragmentCount = totalFragmentCount;

        chromothripsisEvents.push_back(event);
    }

    return chromothripsisEvents;
}

} // namespace sddparser








#endif
