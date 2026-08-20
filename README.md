# SDDparser
Writing a C++ package capable of parsing the Standard for DNA Damage file header and data fields. 

This package is now capable of summarizing the SDD header fields and the data fields into a single summary file. It can also optionally plot a 
Karyogram of the double-strand break locations onto the chromosomes the user passes in the SDD file header.

HOW TO USE:
1. To compile the SDDparser program, navigate to the directory where the Makefile is stored (cd /path/to/SDDparserDirectory/)

2. Ensure g++ is installed on your computer by checking "g++ --version". The terminal should return the version number and license.

3. Ensure you have installed the libcairo2-dev library for plotting the associated Karyogram to the SDD file.

3. Simply type 'make' into your terminal to compile the program 'SDDparser'.

4. To run the SDDparser, type in the command line './SDDparser /path/to/SDDinputFile.txt', where SDDinputFile.txt is the SDD file you want to parse. 

5. The SDD file summary will be stored in a file labeled 'SDDinputFile_summary.txt' in the same directory as the SDDinputFile.

6. If you decide to plot the Karyogram of the SDD file illustrating the locations of the double-strand breaks on each of the chromosomes, 
type in the command-line './SDDparser /path/to/SDDinputFile.txt --karyogram human|other'. '--karyogram' indicates you want to draw
the karyogram of the associated damages, and you must specify either 'human' for human genome centromere positions or 'other' for generic centromere
locations. The output .png file will be stored in the same directory as the SDD input file.

7. IMPORTANT: So far, the karyogram plotter assumes the user will pass the chromosome sizes in the following way for human chromosomes: 1, 2, 3,..., 22, 1,
 2, 3, ..., 22, Y, X. The karyogram plotter also works for more than 46 chromosomes assuming that this chromosome listing practice is upheld. The user must 
specify the homologous chromosomes as well. If the user desires, the cell cycle phase can also be specified in the SDD header, and the karyogram plotter 
will account for if the cell cycle phase is pre-replication (G0 or G1 phase) or post-replication (S, G2, M phase). If the user passes '0' (unspecified) as 
the cell cycle phase, the plotter will assume the cells are in a pre-replication phase. In addition to the chromosome number and sizes passed in the 
SDD header, the karyogram plotter also requires SDD data fields 3, 4, and 6, otherwise the chromosome sizes and damage information will not be present
to draw the karyogram.

WHAT IS SUMMARIZED IN THE SDD FILE SUMMARY:

All important header fields are interpreted and summarized in a text format at the beginning of the summary file.
The header summary is separated into 3 main subsections:
1. Incident Radiation Information - describes the type of radiation particle simulated.
2. Radiation Target Information - describes the cell/nucleus model that was irradiated in the simulation.
3. DNA Damage Information - describes the damage definition for the appearance of Double Strand Breaks, and other exposure information.

The data block is summarized under the subsection 'Chromosome Damages', where the following information is summarized:
1. The number of exposures contained in the SDD file (number of cells irradiated).
2. The number of damage entries (data rows) per exposure.
3. The associated damages per chromosome per exposure (i.e. number base damages, single-strand breaks, and double-strand breaks)
4. Total number of base damages, single-strand breaks, and double-strand breaks over all chromosomes per exposure.
5. The total number of base damages, single-strand breaks, and double-strand breaks over all chromosomes over all exposures.

WHAT IS PLOTTED IN THE SDD FILE KARYOGRAM:
1. All chromosomes passed by the user are plotted with distinct colors for visual clarity. If the user passes in the SDD file a genome that is
not of human origin (i.e. not 46 chromosomes), the user must specify '--karyogram other'.
2. Centromeric positions using the '--karyogram human' option are adjusted based on the chromosome being drawn, whereas non-human genomes will
be drawn using the generic centromere locations at approximately 35% of the chromosomes length. 
3. Double-strand breaks are denoted using black lines spanning slightly more than the width of the chromosome to disitnguish from single-strand breaks.
4. Single-strand breaks are denoted using white lines on the chromosomes and the line lengths are scaled based on the number of single-strand breaks 
within a given damage site (i.e. a given SDD data row) (usually between 0-5 SSBs). 
5. Legend at the top to describe what each symbol signifies in the karyogram.

FURTHER KARYOGRAM PLOTTING WORK:
1. Add a zoom in and out and a panning option to conserve image quality. 
2. Be able to read a Standard for DNA repair (SDR) file and plot the structural variation mutations on the Karyogram.
... and much more.
