#pragma once

#include <cairo/cairo.h>

#include <string>
#include <vector>

class Karyogram
{
public:

    bool generateKaryogram(
        const std::vector<double>& chromosomeSizes,
        const std::string& outputFilename
    );

private:

    void drawChromosome(
        cairo_t* cr,
        double x,
        double y,
        double height
    );
};
