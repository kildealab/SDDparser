#include <cmath>
#include <algorithm>
#include <map>
#include <string>
#include <set>
#include <iostream>

#include "Karyogram.h"
#include "SDRutilities.h"



namespace sddparser
{

// ------------------------------------------------------- //
// Helper functions to determine shape of mutations and
// how they should be drawn
// ------------------------------------------------------- //


// Recognizes the DONOR side of a deletion-translocation event: two
// records in one slot, a "recombined" record (fragments from the
// home strand PLUS exactly one foreign fragment from elsewhere) and
// a "excised" record (a single home-strand-only fragment). Unlike
// isDeletionShape(), the recombined record is allowed foreign
// material, since the swapped-in piece isn't a gap to fill in, it's
// already correctly part of the recombined molecule as drawn.
bool Karyogram::isDeletionTranslocationDonorShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID)
{
    if (records.size() != 2)
    {
        return false;
    }

    const SDRdataRecord* recombinedRecord = nullptr;
    const SDRdataRecord* excisedRecord = nullptr;

    for (const SDRdataRecord* record : records)
    {
        if (!record->linear)
        {
            return false;
        }

        int homeFragmentCount = 0;
        int foreignFragmentCount = 0;

        for (const SDRfragment& fragment : record->fragments)
        {
            if (fragment.oldStrandID == homeOldStrandID)
            {
                if (isReversedFragment(fragment))
                {
                    return false; // The home fragment itself should not be reversed, a real deletion/translocation donor's own material stays in orientation.
                }

                ++homeFragmentCount;
            }
            else
            {
                // Foreign fragment MAY be reversed, represents an
                // inverted exchange (e.g. forming a dicentric/acentric
                // pair), not a reason to reject this shape.
                ++foreignFragmentCount;
            }
        }

        if (homeFragmentCount >= 1 && foreignFragmentCount >= 1 && recombinedRecord == nullptr)
        {
            recombinedRecord = record;
        }
        else if (homeFragmentCount == 1 && foreignFragmentCount == 0 && excisedRecord == nullptr)
        {
            excisedRecord = record;
        }
        else
        {
            return false;
        }
    }

    return recombinedRecord != nullptr && excisedRecord != nullptr;
}















// Recognizes the DONOR side of a deletion-insertion event when it's
// drawn alone in its slot: a single home-strand-only record with 2+
// fragments and at least one real gap between them. The material
// that left this strand shows up as a foreign fragment inside a
// completely different slot's record (the recipient), so there's no
// matching excised record here to pair with, this shape has to be
// recognized from a single record alone.
bool Karyogram::isLoneGapShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID)
{
    if (records.size() != 1)
    {
        return false;
    }

    const SDRdataRecord* record = records[0];

    if (!record->linear)
    {
        return false;
    }

    if (record->fragments.size() < 2)
    {
        return false;
    }

    for (const SDRfragment& fragment : record->fragments)
    {
        if (fragment.oldStrandID != homeOldStrandID)
        {
            return false;
        }

        if (isReversedFragment(fragment))
        {
            return false; // Not this shape - an inversion, not a deletion-insertion donor.
        }
    }

    std::vector<SDRfragment> sortedFragments = record->fragments;

    std::sort(sortedFragments.begin(), sortedFragments.end(), [](const SDRfragment& a, const SDRfragment& b)
    {
        return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
    });

    const double tolerance = 0.001;
    bool hasGap = false;

    for (std::size_t i = 0; i + 1 < sortedFragments.size(); ++i)
    {
        const double higherPosStrandA = std::max(sortedFragments[i].oldStartPosition, sortedFragments[i].oldEndPosition);
        const double lowerPosStrandB = std::min(sortedFragments[i + 1].oldStartPosition, sortedFragments[i + 1].oldEndPosition);

        if (lowerPosStrandB < higherPosStrandA && !approxEqual(higherPosStrandA, lowerPosStrandB, tolerance))
        {
            return false; // Genuine overlap - not a valid shape.
        }

        if (lowerPosStrandB > higherPosStrandA && !approxEqual(higherPosStrandA, lowerPosStrandB, tolerance))
        {
            hasGap = true;
        }
    }

    return hasGap;
}









// Same pairing pattern as isDeletionShape()/isECDNAshape(), but for a
// deletion-inversion: one record has THREE fragments (2 normal + 1
// reversed) from the home strand, with a real gap between two of
// them; a second record supplies the single fragment filling that
// gap. All fragments (both records) must be linear.
bool Karyogram::isDeletionInversionShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID)
{
    if (records.size() != 2)
    {
        return false;
    }

    const SDRdataRecord* threeFragmentRecord = nullptr;
    const SDRdataRecord* oneFragmentRecord = nullptr;

    for (const SDRdataRecord* record : records)
    {
        if (!record->linear)
        {
            return false;
        }

        for (const SDRfragment& fragment : record->fragments)
        {
            if (fragment.oldStrandID != homeOldStrandID)
            {
                return false;
            }
        }

        if (record->fragments.size() == 3 && threeFragmentRecord == nullptr)
        {
            threeFragmentRecord = record;
        }
        else if (record->fragments.size() == 1 && oneFragmentRecord == nullptr)
        {
            oneFragmentRecord = record;
        }
        else
        {
            return false;
        }
    }

    if (threeFragmentRecord == nullptr || oneFragmentRecord == nullptr)
    {
        return false;
    }



    int reversedCount = 0;

    for (const SDRfragment& fragment : threeFragmentRecord->fragments)
    {
        if (fragment.oldStartPosition > fragment.oldEndPosition)
        {
            ++reversedCount;
        }
    }

    if (reversedCount != 1)
    {
        return false;
    }


    std::vector<SDRfragment> sortedThreeFragmentRecord = threeFragmentRecord->fragments;

    std::sort(sortedThreeFragmentRecord.begin(), sortedThreeFragmentRecord.end(), [](const SDRfragment& a, const SDRfragment& b)
    {
        return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
    });

    const SDRfragment& excisedFragment = oneFragmentRecord->fragments[0];

    if (excisedFragment.oldStartPosition > excisedFragment.oldEndPosition)
    {
        return false; // The excised piece itself should not be reversed.
    }

    const double delTolerance = 0.001;
    int gapCount = 0;
    bool gapMatchesExcised = false;

    for (std::size_t i = 0; i + 1 < sortedThreeFragmentRecord.size(); ++i)
    {
        const double higherPosStrandA = std::max(sortedThreeFragmentRecord[i].oldStartPosition, sortedThreeFragmentRecord[i].oldEndPosition);
        const double lowerPosStrandB = std::min(sortedThreeFragmentRecord[i + 1].oldStartPosition, sortedThreeFragmentRecord[i + 1].oldEndPosition);

        if (approxEqual(higherPosStrandA, lowerPosStrandB, delTolerance))
        {
            continue;
        }

        if (lowerPosStrandB <= higherPosStrandA)
        {
            return false; // Overlapping, not a gap - invalid shape.
        }

        ++gapCount;

        if (approxEqual(excisedFragment.oldStartPosition, higherPosStrandA, delTolerance) &&
            approxEqual(excisedFragment.oldEndPosition, lowerPosStrandB, delTolerance))
        {
            gapMatchesExcised = true;
        }
    }

    return gapCount == 1 && gapMatchesExcised;
}















// Identical shape to isDeletionShape(), except the excised (single fragment)
// records must be CIRCULAR rather than linear, and must be acentric -
// an ecDNA fragment carrying a centromere isn't ecDNA. Generalized to N+1
// ecDNA fragments on the same original strand. So instead of requiring exactly
// one fragment per excised record, this checks that the TOTAL
// fragment count across all excised records matches the total number
// of gaps in the remaining record.

bool Karyogram::isECDNAshape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, const SDRmasterHeader& masterHeader)
{
    if (records.size() < 2)
    {
        return false;
    }

    const SDRdataRecord* multipleFragmentRecord = nullptr;
    std::vector<const SDRdataRecord*> excisedRecords;

    for (const SDRdataRecord* record : records)
    {
        for (const SDRfragment& fragment : record->fragments)
        {
            if (fragment.oldStrandID != homeOldStrandID)
            {
                return false;
            }
        }

        if (record->fragments.size() >= 2 && record->linear && multipleFragmentRecord == nullptr)
        {
            multipleFragmentRecord = record;
        }
        else if (!record->linear && !record->fragments.empty())
        {
            bool hasCentromereFlag = false;

            for (const SDRfragment& fragment : record->fragments)
            {
                if (fragment.hasCentromere)
                {
                    hasCentromereFlag = true;
                    break;
                }
            }

            if (hasCentromereFlag)
            {
                return false;
            }

            excisedRecords.push_back(record);
        }
        else
        {
            return false;
        }
    }

    if (multipleFragmentRecord == nullptr || excisedRecords.empty())
    {
        return false;
    }


    std::size_t totalExcisedFragments = 0;

    for (const SDRdataRecord* record : excisedRecords)
    {
        totalExcisedFragments += record->fragments.size();
    }



    std::vector<SDRfragment> sortedFragments = multipleFragmentRecord->fragments;

    std::sort(sortedFragments.begin(), sortedFragments.end(), [](const SDRfragment& a, const SDRfragment& b)
    {
        return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
    });


    const double delTolerance = 0.001;
    std::vector<std::pair<double, double>> gaps;

    const std::size_t sizeIndex = static_cast<std::size_t>(homeOldStrandID);
    const double chromosomeFullSizeMbp = (sizeIndex < masterHeader.intactChromosomeSizes.size()) ? masterHeader.intactChromosomeSizes[sizeIndex] : 0.0;

    // ------------------------------------------------------------- //
    // Determine where the deleted segment came from on the chromosome //
    // ------------------------------------------------------------- //

    // Check if fragment came from the start of the chromosome
    if (chromosomeFullSizeMbp > 0.0)
    {
        const double firstFragmentStart = std::min(sortedFragments.front().oldStartPosition, sortedFragments.front().oldEndPosition);

        if (!approxEqual(firstFragmentStart, 0.0, delTolerance) && firstFragmentStart > 0.0)
        {
            gaps.push_back({0.0, firstFragmentStart});
        }
    }

    // Check if fragment came from the middle of the chromosome
    for (std::size_t i = 0; i + 1 < sortedFragments.size(); i++)
    {
        const double higherPosStrandA = std::max(sortedFragments[i].oldStartPosition, sortedFragments[i].oldEndPosition);
        const double lowerPosStrandB = std::min(sortedFragments[i + 1].oldStartPosition, sortedFragments[i + 1].oldEndPosition);

        if (approxEqual(higherPosStrandA, lowerPosStrandB, delTolerance))
        {
            continue;
        }

        if (lowerPosStrandB <= higherPosStrandA)
        {
            return false; // Overlapping, not a gap - invalid shape.
        }

        gaps.push_back({higherPosStrandA, lowerPosStrandB});
    }

    // Check if fragment came from the end of the chromosome

    if (chromosomeFullSizeMbp > 0.0)
    {
        const double lastFragmentEnd = std::max(sortedFragments.back().oldStartPosition, sortedFragments.back().oldEndPosition);

        if (chromosomeFullSizeMbp > lastFragmentEnd && !approxEqual(lastFragmentEnd, chromosomeFullSizeMbp, delTolerance))
        {
            gaps.push_back({lastFragmentEnd, chromosomeFullSizeMbp});
        }
    }



    if (gaps.size() != totalExcisedFragments)
    {
        return false;
    }


    std::vector<bool> gapClaimed(gaps.size(), false);

    for (const SDRdataRecord* excisedRecord : excisedRecords)
    {
        for (const SDRfragment& excisedFragment : excisedRecord->fragments)
        {
            if (excisedFragment.oldStartPosition > excisedFragment.oldEndPosition)
            {
                return false;
            }


            bool matched = false;

            for (std::size_t i = 0; i < gaps.size(); ++i)
            {
                if (gapClaimed[i])
                {
                    continue;
                }

                if (approxEqual(excisedFragment.oldStartPosition, gaps[i].first, delTolerance) &&
                    approxEqual(excisedFragment.oldEndPosition, gaps[i].second, delTolerance))
                {
                    gapClaimed[i] = true;
                    matched = true;
                    break;
                }
            }

            if (!matched)
            {
                return false;
            }
        }
    }

    return true;

}














// Generalized deletion shape to account for multiple deletions on the same strand.
// Shape looks like this, one record contains N + 1 fragments, where N is the number of
// deletions on that original strand. All N + 1 fragments have the same old strand ID.
// Then this record is followed by N records representing the N deletions, each with one
// fragment in the record. All deletion records have one fragment with the same old strand ID.
bool Karyogram::isDeletionShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, const SDRmasterHeader& masterHeader)
{
    if (records.size() < 2)
    {
        return false;
    }

    const SDRdataRecord* multipleFragmentRecord = nullptr;
    std::vector<const SDRdataRecord*> singleFragmentRecords;

    for (const SDRdataRecord* record : records)
    {

        if(!record->linear)
        {
            return false;
        }

        for (const SDRfragment& fragment : record->fragments)
        {
            if (fragment.oldStrandID != homeOldStrandID)
            {
                return false;
            }
        }

        if (record->fragments.size() >= 2 && multipleFragmentRecord == nullptr)
        {
            multipleFragmentRecord = record;
        }
        else if (record->fragments.size() == 1)
        {
            singleFragmentRecords.push_back(record);
        }
        else
        {
            return false;
        }
    }


    // ------------------------------------------------------------ //
    // Determine if deletion is in the middle or at the ends of a chromosome
    // ------------------------------------------------------------ //
    std::vector<SDRfragment> sortedFragments = multipleFragmentRecord->fragments;

    std::sort(sortedFragments.begin(), sortedFragments.end(), [](const SDRfragment& a, const SDRfragment& b)
    {
        return std::min(a.oldStartPosition, a.oldEndPosition) < std::min(b.oldStartPosition, b.oldEndPosition);
    });

    const double delTolerance = 0.001;
    std::vector<std::pair<double, double>> gaps;


    // A TERMINAL deletion (touching position 0 or the chromosome's true
    // end) leaves no extra "remaining" fragment on that side, so the
    // old rigid N+1 formula only held when every deletion was strictly
    // internal. Checking against the chromosome's true declared length
    // (rather than just gaps between listed fragments) generalizes
    // this to terminal deletions too.
    const std::size_t sizeIndex = static_cast<std::size_t>(homeOldStrandID);
    const double chromosomeFullSizeMbp = (sizeIndex < masterHeader.intactChromosomeSizes.size())
        ? masterHeader.intactChromosomeSizes[sizeIndex]
        : 0.0;

    // Check if deletion is at the start of the chromosome
    if (chromosomeFullSizeMbp > 0.0)
    {
        const double firstFragmentStart = std::min(sortedFragments.front().oldStartPosition, sortedFragments.front().oldEndPosition);

        if (!approxEqual(firstFragmentStart, 0.0, delTolerance) && firstFragmentStart > 0.0)
        {
            gaps.push_back({0.0, firstFragmentStart});
        }
    }



    // Check for deletions in central segments
    for (std::size_t i = 0; i + 1 < sortedFragments.size(); ++i)
    {
        const double higherPosStrandA = std::max(sortedFragments[i].oldStartPosition, sortedFragments[i].oldEndPosition);
        const double lowerPosStrandB = std::min(sortedFragments[i + 1].oldStartPosition, sortedFragments[i + 1].oldEndPosition);

        if (approxEqual(higherPosStrandA, lowerPosStrandB, delTolerance))
        {
            continue;
        }

        if (lowerPosStrandB <= higherPosStrandA)
        {
            return false; // Overlapping, not a gap - invalid shape.
        }

        gaps.push_back({higherPosStrandA, lowerPosStrandB});
    }


    // Check if deletion at the end of the chromosome
    if (chromosomeFullSizeMbp > 0.0)
    {
        const double lastFragmentEnd = std::max(sortedFragments.back().oldStartPosition, sortedFragments.back().oldEndPosition);

        if (chromosomeFullSizeMbp > lastFragmentEnd && !approxEqual(lastFragmentEnd, chromosomeFullSizeMbp, delTolerance))
        {
            gaps.push_back({lastFragmentEnd, chromosomeFullSizeMbp});
        }
    }




    if (gaps.size() != singleFragmentRecords.size())
    {
        return false;
    }

    std::vector<bool> gapClaimed(gaps.size(), false);

    for (const SDRdataRecord* excisedRecord : singleFragmentRecords)
    {
        const SDRfragment& excisedFragment = excisedRecord->fragments[0];

        if (excisedFragment.oldStartPosition > excisedFragment.oldEndPosition)
        {
            return false; // Excised pieces should not be reversed.
        }

        bool matched = false;

        for (std::size_t i = 0; i < gaps.size(); ++i)
        {
            if (gapClaimed[i])
            {
                continue;
            }

            if (approxEqual(excisedFragment.oldStartPosition, gaps[i].first, delTolerance) &&
                approxEqual(excisedFragment.oldEndPosition, gaps[i].second, delTolerance))
            {
                gapClaimed[i] = true;
                matched = true;
                break;
            }
        }

        if (!matched)
        {
            return false;
        }
    }

    return true;
}

} // namespace sddparser
