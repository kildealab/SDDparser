# SDDparser
Writing a C++ package capable of parsing the Standard for DNA Damage file header and data fields. 

This package is now capable of summarizing the SDD header fields and the data fields into a single summary file. It can also optionally plot a 
Karyogram of the double-strand break and single-strand break locations onto the chromosomes the user passes in the SDD file header. Current work
is being done to accomodate Standard for DNA Repair (SDR) file parsing and karyogram plotting.

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

7. IMPORTANT: The karyogram can now handle the user passing 'Chromosome sizes' in the SDD header in the following three ways (example for human chromosomes) :

	a. Split homolog chromosome sizes layout: 1,2,3,...,22,1,2,3,...,22,Y,X.

	b. Adjacent homolog chromosome sizes layout: 1,1,2,2,3,3,...,22,22,Y,X.

	c. Non-homologous/haploid chromosome sizes layout: 1,2,3,...,22,Y,X.

In the SDD chromosome sizes header, please specify Y chromosome size before X. If the user passes two X chromosomes, the second one will be labeled as Y
and will have an incorrect centromere position in the karyogram. 
The plotter also works for more than 46 chromosomes, but the '--karyogram other' option must be specified by the user.  
If the user desires, the cell cycle phase can also be specified in the SDD header, and the karyogram plotter will account for if the cell cycle phase 
is pre-replication (G0 or G1 phase) or post-replication (S, G2, M phase). If the user passes '0' (unspecified) as the cell cycle phase, the plotter 
will assume the cells are in a pre-replication phase. In addition to the chromosome number and sizes passed in the SDD header, the karyogram plotter 
also requires SDD data fields 3, 4, and 6, otherwise the chromosome sizes and damage information will not be present to draw the karyogram. 

****PLEASE ENSURE THE CHROMOSOME IDS IN DATA FIELD 3 CORRESPOND TO THE CORRECT CHROMOSOME SIZE INDEX LISTED IN THE SDD HEADER, 
OTHERWISE THE KARYOGRAM PLOTTER WILL BE INCORRECT****

8. Functionality has been added to parse an SDR file using the following command './SDDparser -sdr ./path/to/SDRinputFile.txt'. The summary of the SDR header
and subheader are returned for each cell, as well as a summary of the number of mutations present in the SDR file (so far long deletions, balanced inversion, 
and balanced translocations). For now, karyogram plotting is not supported for depicting genomic rearrangements from DNA structural variations, but this
should be implemented soon. 

To check if everything works correctly, run the following command:
'./SDDparser exampleSDD.txt --karyogram human'

The output to the terminal should be:
'SDD header parsed successfully.
Summary written to: "exampleSDD_summary.txt"
Karyogram generated successfully: "exampleSDD_karyogram_human.png"'
The exampleSDD_summary.txt file and exampleSDD_karyogram_human.png files
should resemble the ones shown in the attached files of this repository.

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
2. Centromeric positions using the '--karyogram human' option are realistic and adjusted based on the chromosome being drawn, whereas 
non-human genomes will be drawn using the generic centromere locations at approximately 35% of the chromosomes length. 
3. Double-strand breaks are denoted using black lines spanning slightly more than the width of the chromosome to disitnguish from single-strand breaks.
4. Single-strand breaks are denoted using white lines on the chromosomes and the line lengths are scaled based on the number of single-strand breaks 
within a given damage site (i.e. a given SDD data row) (usually between 0-5 SSBs). 
5. Summary at the top to describe the dose or fluence in the SDD file, the cell cycle phase, and the number of single- and double-strand breaks.
6. Legend at the bottom to describe what each symbol signifies in the karyogram.

WHAT IS SUMMARIZED IN THE SDR FILE SUMMARY:
1. SDR version, author, associated SDD file that produced the SDR file from the MEDRAS-MC output.
2. Number of chromosomes listed and their respective sizes in mega base pairs.
3. Per Cell summary of the cell subheader (number of double-strand breaks and misrepairs).
4. Number of mutations per each type (so far only long deletions, balanced inversions, and balanced translocations are supported). More support
for more complex structural variations will be implemented. 

FURTHER KARYOGRAM PLOTTING WORK:
1. Add a zoom in and out and a panning option to conserve image quality. 
2. Be able to read a Standard for DNA repair (SDR) file and plot the structural variation mutations on the Karyogram.
3. Add functionality to be able to receive two X chromosomes instead of a Y and X chromosome and be able to plot them on the karyogram. 
4. Depict the increase in number of single-strand and double-strand breaks as a function of dose using the Karyogram for SDD files with different doses listed,
... and much more.
