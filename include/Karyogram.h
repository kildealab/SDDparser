#pragma once

#include <cairo/cairo.h>

#include <string>
#include <vector>


struct RGB
{
    double r;
    double g;
    double b;
};

class Karyogram
{
public:

    bool generateKaryogram(
        const std::vector<double>& chromosomeSizes,
	const std::vector<double>& cellCyclePhase, // Check to know if Karyogram should draw single or double chromatid arms.
        const std::string& outputFilename
    );

private:

    RGB generateChromosomeColor(int chromosomeNumber); 
		

    void drawChromosome(
        cairo_t* cr,
        double x,
        double y,
        double height,
	RGB color
    );
};
