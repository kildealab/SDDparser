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
	const SDRmasterHeader& masterHeader,		// One master header per SDR file containing most importantly the sizes of the original chromosomes before mutations
	const SDRsubHeader& subHeader,			// One subheader per cell in the SDR file, iterate through each subheader for each cell
	bool humanGenome,				// Check if user specified human genome to adjust centromere positions on karyogram
	const std::string& outputFilename		// Output file names for the summary file and the karyogram drawing.
    );


private:

    RGB generateChromosomeColor(int chromosomeNumber, int totalChromosomes); 	// Function to generate chromosome colors depending on the number of chromosomes, and to choose colors that are different enough between successive chromosomes. 
    RGB getColorForOriginalStrand(int oldStrandID, const SDRmasterHeader& masterHeader);			// For SDR tracking of chromosome colors during rearrangements

    // Helper functions to build the mutated segments onto the original chromosome strands, adjusting the color depending on the origin of the mutation, or white for deletions
    std::vector<PaintedSegment> buildPaintedSegments(const SDRdataRecord& record, bool humanGenome, const SDRmasterHeader& masterHeader);
    std::vector<PaintedSegment> buildDeletionRemainingSegments(const SDRdataRecord& record, int homeOldStrandID, bool humanGenome, const SDRmasterHeader& masterHeader);
    std::vector<PaintedSegment> buildMixedStrandSegmentsWithGaps(const SDRdataRecord& record, int homeOldStrandID, bool humanGenome, const SDRmasterHeader& masterHeader, double& outTotalLengthMbp);

    void drawChromosome(cairo_t* cr, double x, double y, double height, double width, RGB color, double centromereStart, double centromereEnd); // The main draw chromosome function, accounting for individual chromosome sizes and centromere ranges and locations.
    void drawPaintedChromosome(cairo_t* cr, double x, double y, double height, double width, const std::vector<PaintedSegment>& segments, bool roundTopCap = true, bool roundBottomCap = true, bool drawOutline = true);
    void drawStackedMutations(cairo_t* cr, const std::vector<const SDRdataRecord*>& records, double slotCenterX, double posY, double chromosomeWidth, double maxLengthMbp, double maxRenderHeight, bool humanGenome, const SDRmasterHeader& masterHeader, int homeOldStrandID);

    std::vector<DamageLocation> getDoubleStrandBreaks(const std::vector<Exposure>& exposures); // Use to obtain the stored double strand break locations in each exposure and determine their coordinates on the Karyogram. 
    std::vector<DamageLocation> getSingleStrandBreaks(const std::vector<Exposure>& exposures);

    double getDamageFraction(const DamageLocation& damage, double chromosomeSize);		// Convert damage locations in base pairs to a fractional length along the chromosome for easy Karyogram pixel conversion.

    void drawDoubleStrandBreakMarker(cairo_t* cr, double x, double y, double chromosomeWidth);				// Functions to draw the damage locations at a given x and y coordinate on the Karyogram, the entire width of the drawn chromosome.
    void drawSingleStrandBreakMarker(cairo_t* cr, double x, double y, double markerLength);
    void drawInversionChevron(cairo_t* cr, double centerX, double centerY, double size);	// Draw symbol for SDR file inversions
    void drawCircularFragment(cairo_t* cr, double centerX, double centerY, double diameter, RGB color, bool drawOutline = true);	// Draw symbol for SDR file ecDNA 

    const CentromerePosition* getHumanCentromere(int chromosomeID);		// Returns a chromosome's corresponding centromere start and end locations to draw the centromere ellipse on the karyogram.
    bool getCentromereForOriginalStrand(int oldStrandID, bool humanGenome, const SDRmasterHeader& masterHeader, double&centromereStartBP, double& centromereEndBP);			// Check if original strand has centromere in the given range

    // Helper functions to modify the chromosome heights and widths to adjust for mutated content being carried over into this slot.
    double computeSDRbarHeight(double lengthMbp, double maxLengthMbp, double maxRenderHeight);
    double computeMaxBarHeight(const std::vector<const SDRdataRecord*>& records, double maxLengthMbp, double maxRenderHeight);
    double computeSlotWidth(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, double chromosomeWidth, double maxLengthMbp, double maxRenderHeight, const SDRmasterHeader& masterHeader, double& outColumn2Width);
    int computeExcisedColumnLayout(const std::vector<const SDRdataRecord*>& excisedRecords, double maxLengthMbp, double maxRenderHeight, double maxColumnHeight, double verticalGap, std::vector<double>& outHeights);

    void drawLegend(cairo_t* cr, double legendY);						// Function to draw the SDD karyogram legend at the bottom.
    void drawSDDsummary(cairo_t* cr, const std::vector<double>& cellCyclePhase, const std::vector<Exposure>& exposures, const std::vector<double>& doseOrFluence, const std::vector<int>& incidentParticles);						// Function to draw the karyogram summary box at the top.
    void drawSDRsummary(cairo_t* cr, const SDRmasterHeader& masterHeader, const SDRsubHeader& subHeader);	// Function to draw the top SDR karyogram summary box summarizing which mutations were in the SDR file
    void drawSDRlegend(cairo_t* cr, double legendY);						// Draw legend at bottom of the SDR karyogram which has different symbols than the SDD file karyogram

    ChromosomeLayout determineChromosomeLayout(const std::vector<double>& chromosomeSizes);	// Function that will determine how the user passed the chromosomes.

    // Function to check if the SDR data entries are depciting mutations, or are the original intact strand information
    std::vector<const SDRdataRecord*> filterBaselineIfMutated(const std::vector<const SDRdataRecord*>& records, int numOriginalStrands);

    // Functions to alter the drawing of certain mutations depending on the specific mutation involved
    bool isDeletionShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID);		// For long deletions
    bool isECDNAshape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID);			// For ecDNA mutations
    bool isDeletionInversionShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID);	// For deletion-inversion mutations
    bool isDeletionTranslocationDonorShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID);	// For deletion-translocation mutations
    bool isLoneGapShape(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID);			// For deletion-insertion mutations

    // Helper functions to stack mutated chromosome segments of different colors on top of the original segments and remove the correct amount of 
    // chromosome if a deletion occurs.
    SDRdataRecord concatenateRecordsForDrawing(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID);
    std::vector<const SDRdataRecord*> clusterRecordsForDrawing(const std::vector<const SDRdataRecord*>& records, int homeOldStrandID, std::vector<SDRdataRecord>& mergedStorage);

    // If intact chromosome strands are not listed in SDR data record, synthesize them using the intact chromosome sizes.
    std::vector<const SDRdataRecord*> synthesizeIntactRecordIfMissing(int oldStrandID, int cellID, const SDRmasterHeader& masterHeader, std::vector<SDRdataRecord>& intactRecordStorage);

};

#endif
