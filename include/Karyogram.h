#ifndef KARYOGRAM_H
#define KARYOGRAM_H

#include <cairo/cairo.h>

#include <string>
#include <vector>
#include <utility>

#include "SDDtypes.h"
#include "SDRtypes.h"

struct RGB                                      // Use to define chromosome colors using Red, Green, and Blue codes.
{
    double r;
    double g;
    double b;
};


struct PaintedSegment
{
    double startFraction;                       // Position along new strand as fraction of total length
    double endFraction;
    RGB color;                                  // Color of original chromosome this segment came from
    bool hasCentromere;                         // Does the centromere exist in this segment
    bool isReversed;                            // Check for inversions to draw chevron inversion marker.
    double centromereStartFraction;
    double centromereEndFraction;
};


enum class ChromosomeLayout                     // Determine how the user passed the chromosome sizes to modify karyogram plotting logic.
{
    NON_HOMOLOGOUS,       // 1,2,3,...,22,Y,X
    ADJACENT_HOMOLOGS,    // 1,1,2,2,...,22,22,Y,X
    SPLIT_HOMOLOGS        // 1,2,...,22,1,2,...,22,Y,X
};


struct CentromerePosition                       // Map the centromere start and end positions to the corresponding chromosome for the human genome alone.
{
    std::string chromosome;                     // 1-22, X, Y
    long long start;                            // Start position in base pairs
    long long end;                              // End position in base pairs

};




// Use Karyogram for all illustrations of chromosome damages
class Karyogram
{
public:

    bool generateKaryogram(				// Function that generates the overall Karyogram structure and calls individual functions to draw each chromosome based on sizes.
        const std::vector<double>& chromosomeSizes,	// Use to scale the chromosome heights and convert to number of pixels on the Karyogram.
	const std::vector<double>& cellCyclePhase, 	// Look at SDD header 'Cell cycle phase' to determine if chromosomes should be replicated or not in the final karyogram.
	const std::vector<double>& doseOrFluence,	// Look at SDD header for dose or fluence to be summarized in Karyogram.
	const std::vector<int>& incidentParticles,	// Look at SDD header for 'Incident Particles' to be summarized in karyogram.
	const std::vector<Exposure>& exposures,		// Loop through Exposures object to label chromosome IDs onto the karyogram
	bool humanGenome,				// Check if the user passed the argument --karyogram human|other. if humanGenome = true, adjust the centromere locations. If humanGenome = false, use generic centromere locations.
        const std::string& outputFilename		// Specify generic output file name combining the input file prefix with the suffix "_karyogram.png"
    );

    bool generateSDRkaryogram(				// Function to generate the karyogram of the rearrangements in an SDR file
	const SDRmasterHeader& masterHeader,
	const SDRsubHeader& subHeader,
	bool humanGenome,
	const std::string& outputFilename
    );

private:

    RGB generateChromosomeColor(int chromosomeNumber, int totalChromosomes); 	// Function to generate chromosome colors depending on the number of chromosomes, and to choose colors that are different enough between successive chromosomes. 
    RGB getColorForOriginalStrand(int oldStrandID, const SDRmasterHeader& masterHeader);			// For SDR tracking of chromosome colors during rearrangements

    std::vector<PaintedSegment> buildPaintedSegments(const SDRdataRecord& record, bool humanGenome, const SDRmasterHeader& masterHeader);

    void drawChromosome(cairo_t* cr, double x, double y, double height, double width, RGB color, double centromereStart, double centromereEnd); // The main draw chromosome function, accounting for individual chromosome sizes and centromere ranges and locations.
    void drawPaintedChromosome(cairo_t* cr, double x, double y, double height, double width, const std::vector<PaintedSegment>& segments);
    void drawStackedMutations(cairo_t* cr, const std::vector<const SDRdataRecord*>& records, double slotCenterX, double posY, double chromosomeWidth, double maxLengthMbp, double maxRenderHeight, bool humanGenome, const SDRmasterHeader& masterHeader);

    std::vector<DamageLocation> getDoubleStrandBreaks(const std::vector<Exposure>& exposures); // Use to obtain the stored double strand break locations in each exposure and determine their coordinates on the Karyogram. 
    std::vector<DamageLocation> getSingleStrandBreaks(const std::vector<Exposure>& exposures);

    double getDamageFraction(const DamageLocation& damage, double chromosomeSize);		// Convert damage locations in base pairs to a fractional length along the chromosome for easy Karyogram pixel conversion.

    void drawDoubleStrandBreakMarker(cairo_t* cr, double x, double y, double chromosomeWidth);				// Function to draw the damage locations at a given x and y coordinate on the Karyogram, the entire width of the drawn chromosome.
    void drawSingleStrandBreakMarker(cairo_t* cr, double x, double y, double markerLength);
    void drawInversionChevron(cairo_t* cr, double centerX, double centerY, double size);

    const CentromerePosition* getHumanCentromere(int chromosomeID);		// Returns a chromosome's corresponding centromere start and end locations to draw the centromere ellipse on the karyogram.
    bool getCentromereForOriginalStrand(int oldStrandID, bool humanGenome, const SDRmasterHeader& masterHeader, double&centromereStartBP, double& centromereEndBP);			// Check if original strand has centromere in the given range

    double computeSDRbarHeight(double lengthMbp, double maxLengthMbp, double maxRenderHeight);
    double computeMaxBarHeight(const std::vector<const SDRdataRecord*>& records, double maxLengthMbp, double maxRenderHeight);

    void drawLegend(cairo_t* cr, double legendY);						// Function to draw the karyogram legend at the bottom.
    void drawSDDsummary(cairo_t* cr, const std::vector<double>& cellCyclePhase, const std::vector<Exposure>& exposures, const std::vector<double>& doseOrFluence, const std::vector<int>& incidentParticles);						// Function to draw the karyogram summary box at the top.

    ChromosomeLayout determineChromosomeLayout(const std::vector<double>& chromosomeSizes);	// Function that will determine how the user passed the chromosomes.

    std::vector<const SDRdataRecord*> filterBaselineIfMutated(const std::vector<const SDRdataRecord*>& records, int numOriginalStrands);

};

#endif
