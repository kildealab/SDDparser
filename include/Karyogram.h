#pragma once

#include <cairo/cairo.h>

#include <string>
#include <vector>
#include <utility>

#include "SDDtypes.h"

struct RGB					// Use to define chromosome colors using Red, Green, and Blue codes.
{
    double r;
    double g;
    double b;
};

struct ChromosomeGeometry			// Store the geometry to draw corresponding chromosome damages at a given x and y coordinate.
{
    int chromosomeNumber = 0;
//    int homologNumber = 0;
    int chromatidNumber = 0;
    double x = 0.0;
    double y = 0.0;
    double height = 0.0;
};

struct CentromerePosition			// Map the centromere start and end positions to the corresponding chromosome for the human genome alone.
{
    std::string chromosome; 			// 1-22, X, Y
    long long start; 				// Start position in base pairs
    long long end; 				// End position in base pairs

};

// Use Karyogram for all illustrations of chromosome damages
class Karyogram
{
public:

    bool generateKaryogram(				// Function that generates the overall Karyogram structure and calls individual functions to draw each chromosome based on sizes.
        const std::vector<double>& chromosomeSizes,	// Use to scale the chromosome heights and convert to number of pixels on the Karyogram.
	const std::vector<double>& cellCyclePhase, 	// Look at SDD header 'Cell cycle phase' to determine if chromosomes should be replicated or not in the final karyogram.
	const std::vector<Exposure>& exposures,		// Loop through Exposures object to label chromosome IDs onto the karyogram
	bool humanGenome,				// Check if the user passed the argument --karyogram human|other. if humanGenome = true, adjust the centromere locations. If humanGenome = false, use generic centromere locations.
        const std::string& outputFilename		// Specify generic output file name combining the input file prefix with the suffix "_karyogram.png"
    );

private:

    RGB generateChromosomeColor(int chromosomeNumber); 	// Function to generate chromosome colors depending on the number of chromosomes, and to choose colors that are different enough between successive chromosomes. 

    void drawChromosome(cairo_t* cr, double x, double y, double height, RGB color, double centromereStart, double centromereEnd); // The main draw chromosome function, accounting for individual chromosome sizes and centromere ranges and locations.

    std::vector<DamageLocation> getDoubleStrandBreaks(const std::vector<Exposure>& exposures); // Use to obtain the stored double strand break locations in each exposure and determine their coordinates on the Karyogram. 

    std::vector<DamageLocation> getSingleStrandBreaks(const std::vector<Exposure>& exposures);

    double getDamageFraction(const DamageLocation& damage, double chromosomeSize);		// Convert damage locations in base pairs to a fractional length along the chromosome for easy Karyogram pixel conversion.

    void drawDoubleStrandBreakMarker(cairo_t* cr, double x, double y, double chromosomeWidth);				// Function to draw the damage locations at a given x and y coordinate on the Karyogram, the entire width of the drawn chromosome.

    void drawSingleStrandBreakMarker(cairo_t* cr, double x, double y, int numSingleStrandBreaks);

    const CentromerePosition* getHumanCentromere(int chromosomeID);		// Returns a chromosome's corresponding centromere start and end locations to draw the centromere ellipse on the karyogram.

    void drawLegend(cairo_t* cr);						// Function to draw the karyogram legend.
};
