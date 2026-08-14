#pragma once

#include <cairo/cairo.h>

#include <string>
#include <vector>
#include <utility>

#include "SDDtypes.h"

struct RGB
{
    double r;
    double g;
    double b;
};

struct ChromosomeGeometry
{
    int chromosomeNumber = 0;
//    int homologNumber = 0;
    int chromatidNumber = 0;

    double x = 0.0;
    double y = 0.0;
    double height = 0.0;
};

class Karyogram
{
public:

    bool generateKaryogram(
        const std::vector<double>& chromosomeSizes,
	const std::vector<double>& cellCyclePhase, // Check to know if Karyogram should draw single or double chromatid arms.
	const std::vector<Exposure>& exposures,
        const std::string& outputFilename
    );

private:

    RGB generateChromosomeColor(int chromosomeNumber); 

    void drawChromosome(cairo_t* cr, double x, double y, double height, RGB color);

    std::vector<DSBlocation>
    getDoubleStrandBreaks(const std::vector<Exposure>& exposures);

    double getDamageFraction(const DSBlocation& dsb, double chromosomeSize);

    void drawDamageMarker(cairo_t* cr, double x, double y);

};
