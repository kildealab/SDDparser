# SDDparser
Writing a C++ package capable of parsing the Standard for DNA Damage file header and data fields. 

This package is now capable of summarizing the SDD header fields and the data fields into a single summary file.

HOW TO USE:
1. To compile the SDDparser program, navigate to the directory where the files are stored (cd /path/to/SDDparserDirectory/)

2. Ensure g++ is installed on your computer by checking "g++ --version". The terminal should return the version number and license.

3. Simply type 'make' into your terminal to compile the program 'SDDparser'.

4. To run the SDDparser, type in the command line './SDDparser /path/to/SDDinputFile.txt', where SDDinputFile.txt is the SDD file you want to parse. 

5. The SDD file summary will be stored in a file labeled 'SDDinputFile_summary.txt' in the same directory as the SDDinputFile.


WHAT IS SUMMARIZED:

All important header fields are interpreted and summarized in a text format at the beginning of the summary file.
The header summary is separated into 3 main subsections:
1. Incident Radiation Information - describes the type of radiation particle simulated.
2. Radiation Target Information - describes the cell/nucleus model that was irradiated in the simulation.
3. DNA Damage Information - describes the damage definition for the appearance of Double Strand Breaks, and other exposure information.

The data block is summarized under the subsection 'Chromosome Damages', where the following information is summarized:
1. The number of exposures contained in the SDD file.
2. The number of damage entries (data rows) per exposure.
3. The associated damages per chromosome per exposure (i.e. number base damages, single-strand breaks, and double-strand breaks)
4. Total number of base damages, single-strand breaks, and double-strand breaks over all chromosomes per exposure.
5. The total number of base damages, single-strand breaks, and double-strand breaks over all chromosomes over all exposures.

Further work is currently being done to produce a visualization of the double strand breaks for each chromosome in a Karyotype diagram.
The visualization will probably be done in C++ using the libcairo library.
