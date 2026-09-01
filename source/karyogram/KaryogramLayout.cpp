#include <cairo/cairo.h>
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>

#include "Karyogram.h"


ChromosomeLayout Karyogram::determineChromosomeLayout(
    const std::vector<double>& chromosomeSizes)							// Access chromosome sizes to determine how they were passed by the user, adjacent-homologs, split-homologs, non-homologous.
{
    if (chromosomeSizes.empty())
    {
        std::cerr << "Require Chromosome sizes in mega base pairs to plot karyogram.\n";
        return ChromosomeLayout::NON_HOMOLOGOUS;
    }

    const int chromosomeCount = static_cast<int>(chromosomeSizes[0]);						// The first entry is the number of chromosomes.

    if (chromosomeCount < 2 || chromosomeSizes.size() != static_cast<size_t>(chromosomeCount + 1))		// Need at least two chromosomes to determine a mapping.
    {
        return ChromosomeLayout::NON_HOMOLOGOUS;
    }

    // --------------------------------------------------
    // Helper for comparing chromosome sizes.
    // Add a tolerance to the homologous chromosome sizes
    // of 1e-6 * chromosome size <= 250 bases.
    // --------------------------------------------------
    const auto approximatelyEqual = [](double a, double b)
        {
            const double tolerance = 1e-6;
            return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
        };

    // --------------------------------------------------
    // Case 1: Adjacent homologs
    //
    // chromosomeSizes:
    // [count, 1, 1, 2, 2, 3, 3, ...]
    // --------------------------------------------------

    if (chromosomeCount >= 2 && approximatelyEqual(chromosomeSizes[1], chromosomeSizes[2]))			// If neighbouring chromosome sizes are approximately equal --> adjacent homologous layout
    {
        return ChromosomeLayout::ADJACENT_HOMOLOGS;
    }


    // --------------------------------------------------
    // Case 2: Split homologs
    //
    // chromosomeSizes:
    // [count,1,2,3,...,N,1,2,3,...,N]
    //
    // Compare the first chromosome with the chromosome
    // halfway through the chromosome list.
    // --------------------------------------------------

    if (chromosomeCount >= 4 && (chromosomeCount - 2) % 2 == 0)
    {
        const int homologousSetSize = (chromosomeCount - 2) / 2;

        const size_t firstHomologIndex = 1;

        const size_t secondHomologIndex = static_cast<size_t>(homologousSetSize + 1);

        if (approximatelyEqual(chromosomeSizes[firstHomologIndex], chromosomeSizes[secondHomologIndex]))	// If the first chromosome size is approx equal to the chromosome size halfway through the list (minus the two sex chromosomes) --> split homologous set
        {
            return ChromosomeLayout::SPLIT_HOMOLOGS;
        }
    }


    // --------------------------------------------------
    // Case 3: Non-homologous
    //
    // chromosomeSizes:
    // [count,1,2,3,4,...,N]
    // --------------------------------------------------

    return ChromosomeLayout::NON_HOMOLOGOUS;									// Default case.

}
















const std::vector<CentromerePosition> humanCentromeres =                   // Vector storing the associated chromosome ID centromere positions to
{                                                                          // draw the centromere location if the human option of '--karyotype' is passed by the user in the command-line
    {"1", 121535434, 124535434},                                           // Each vector element contains {'Chromosome ID', Centromere Start Position (bp), Centromere End Position (bp)}
    {"2",  92326171,  95326171},
    {"3",  90504854,  93504854},
    {"4",  49660117,  52660117},
    {"5",  46405641,  49405641},
    {"6",  58830166,  61830166},
    {"7",  58054331,  61054331},
    {"8",  43838887,  46838887},
    {"9",  47367679,  50367679},
    {"10", 39254935,  42245935},
    {"11", 51644205,  54644205},
    {"12", 34856694,  37856694},
    {"13", 16000000,  19000000},
    {"14", 16000000,  19000000},
    {"15", 17000000,  20000000},
    {"16", 35335801,  38335801},
    {"17", 22263006,  25263006},
    {"18", 15460898,  18460898},
    {"19", 24681782,  27681782},
    {"20", 26369569,  29369569},
    {"21", 11288129,  14288129},
    {"22", 13000000,  16000000},
    {"Y", 10100001, 10327000},
    {"X", 58632012, 61632012}
};










const CentromerePosition* Karyogram::getHumanCentromere(				// Function to return the centromere start and end positions for each chromosome number. 
    int chromosomeID)									// Obtain centromere ranges for a given chromosome ID

{
    std::string chromosome;

    if (chromosomeID >= 1 && chromosomeID <= 22)
    {
        chromosome = std::to_string(chromosomeID);
    }
    else if (chromosomeID >= 23 && chromosomeID <= 44)						// Chromosome number 23-44 correspond to the homologs of chromosomes 1-22. Y is chromosome 45 and X is chromosome 46 --> MUST MAINTAIN THIS ORDERING
    {
        // Second homologous set.
        chromosome = std::to_string(chromosomeID - 22);
    }
    else if (chromosomeID == 45)                                           			// Chromosome Number 45 = Y chromosome
    {
        chromosome = "Y";
    }
    else if (chromosomeID == 46)                                           			// Chromosome Number 46 = X chromosome
    {
        chromosome = "X";
    }
    else
    {
        return nullptr;                                                    			// If other Chromosome Numbers encountered, return nullptr.
    }

    for (const auto& centromere : humanCentromeres)                        			// Return the centromere start and end locations for each chromosome ID.
    {
        if (centromere.chromosome == chromosome)
        {
            return &centromere;
        }
    }

    return nullptr;
}











// For SDR file, looks up centromere location (in bp) for the original strand ID, return false if no centromere information
bool Karyogram::getCentromereForOriginalStrand(
    int oldStrandID,										// Need old strand ID to associate a given centromere range
    bool humanGenome,										// Did user pass a human genome? If not, use generic centromere regions at 50% of the chromosome length
    const SDRmasterHeader& masterHeader,							// Access SDR master header for chromosome sizes to decide where to draw the centromeres
    double& centromereStartBP,									// Access centromere range start in BP, converted later to pixels
    double& centromereEndBP)									// Access centromere range end in BP, converted later to pixels
{

    if (humanGenome)
    {
        // SDR has 1-indexing for strand IDs, same as getHumanCentromere
        const CentromerePosition* centromere = getHumanCentromere(oldStrandID);

        if (centromere == nullptr)
        {
            return false;
        }

        centromereStartBP = static_cast<double>(centromere->start);
        centromereEndBP = static_cast<double>(centromere->end);

        return true;
    }
    else
    {
        const std::size_t sizeIndex = static_cast<std::size_t>(oldStrandID); // 0th index is number of chromosomes listed

        if (sizeIndex >= masterHeader.intactChromosomeSizes.size())
        {
            return false;
        }

        const double chromosomeSizeBP = masterHeader.intactChromosomeSizes[sizeIndex] * 1000000.0;  // Convert chromosome size in Mbp to bp

	// Generic centromere start and end locations halfway through the chromosome length for non-human chromosomes.
        centromereStartBP = 0.49 * chromosomeSizeBP;
        centromereEndBP = 0.51 * chromosomeSizeBP;

        return true;
    }
}


