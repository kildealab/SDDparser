# SDDparser
Writing a C++ package capable of parsing the Standard for DNA Damage (SDD) file header and data fields,
and the Standard for DNA Repair (SDR) header and data fields, plotting the associated damages and mutations on a karyogram. 

This package is now capable of summarizing the SDD header fields and the data fields into a single summary file. It can also optionally plot a 
karyogram of the double-strand break and single-strand break locations onto the chromosomes the user passes in the SDD file header. 
This package can also accomodate the SDR file format in the same way as the SDD file, by summarizing the mutations specified
in the SDR file for each cell, as well as optionally plotting the karyogram of the mutations specified in the SDR file.

HOW TO USE:
1. To compile the SDDparser program, navigate to the directory where the Makefile is stored: 'cd /path/to/SDDparserDirectory/'.

2. Ensure g++ is installed on your computer by checking "g++ --version". The terminal should return the version number and license.

3. Ensure you have installed the libcairo2-dev library for plotting the associated Karyogram to the SDD file.

3. Simply type 'make' into your terminal to compile the program 'SDDparser'. No warnings or errors should appear. 

4. To run the SDDparser and generate a summary, type in the command line './SDDparser /path/to/SDDinputFile.txt', where SDDinputFile.txt 
is the SDD file you want to parse. 

5. The SDD file summary will be stored in a file labeled 'SDDinputFile_summary.txt' in the same directory as the SDDinputFile.

6. If you decide to plot the Karyogram of the SDD file illustrating the locations of the single- and double-strand breaks on each of the chromosomes, 
type in the command-line './SDDparser /path/to/SDDinputFile.txt --karyogram human|other'. '--karyogram' indicates you want to draw
the karyogram of the associated damages, and you must specify either 'human' for human genome centromere positions or 'other' for generic centromere
locations. The output .png file will be stored in the same directory as the SDD input file.

7. IMPORTANT: The karyogram can now handle the user passing 'Chromosome sizes' in the SDD and SDR header in the following three ways 
(example for human chromosomes) :

	a. Split homolog chromosome sizes layout: 1,2,3,...,22,1,2,3,...,22,Y,X.

	b. Adjacent homolog chromosome sizes layout: 1,1,2,2,3,3,...,22,22,Y,X.

	c. Non-homologous/haploid chromosome sizes layout: 1,2,3,...,22,Y,X.

In the SDD chromosome sizes header, please specify Y chromosome size before X. If the user passes two X chromosomes, the second one will be labeled as Y
and will have an incorrect centromere position in the karyogram. 
The plotter also works for more than 46 chromosomes, but the '--karyogram other' option must be specified by the user if that is the case.
If the user desires, the cell cycle phase can also be specified in the SDD header, and the karyogram plotter will account for if the cell cycle phase 
is pre-replication (G0 or G1 phase) or post-replication (S, G2, M phase). If the user passes '0' (unspecified) as the cell cycle phase, the plotter 
will assume the cells are in a pre-replication phase. In addition to the chromosome number and sizes passed in the SDD header, the karyogram plotter 
also requires SDD data fields 3, 4, and 6, otherwise the chromosome sizes and damage information will not be present to draw the karyogram. 

****PLEASE ENSURE THE CHROMOSOME IDS IN DATA FIELD 3 CORRESPOND TO THE CORRECT CHROMOSOME SIZE INDEX LISTED IN THE SDD HEADER, 
OTHERWISE THE KARYOGRAM PLOTTER WILL BE INCORRECT****

8. Functionality has been added to parse an SDR file using the following command './SDDparser -sdr ./path/to/SDRinputFile.txt'. The summary of the SDR 
header and subheader are returned for each cell, as well as a summary of the number of mutations present in the SDR file. DO NOT FORGET THE '-sdr' after
./SDDparser, otherwise the program will assume you are feeding an SDD file as input, and exit with an error. 

9. The user can also now optionally plot a karyogram of the mutations described in an SDR file using the following command:
'/path/to/SDDparser -sdr path/to/SDRinputFile.txt --karyogram human|other'. 

To check if the SDD and SDR parsing, summary, and karyogram plotting work correctly, run the following commands:

a) './SDDparser exampleSDD.txt --karyogram human'

The output to the terminal should be:
'SDD header parsed successfully.
Summary written to: "exampleSDD_summary.txt"
Karyogram generated successfully: "exampleSDD_karyogram_human.png"'

(The exampleSDD_summary.txt file and exampleSDD_karyogram_human.png files should resemble the ones shown in the attached files of this repository in the
examples/ folder).

b) './SDDparser -sdr exampleSDR.txt --karyogram human'

The output to the terminal should be:
'Summary written to: "exampleSDR_summary.txt"
Karyogram generated successfully: "exampleSDR_cell0_karyogram_human.png"'

The exampleSDR_summary.txt file and exampleSDR_cell0_karyogram_human.png files should resemble the ones shown in the attached files of this repository in
the examples/ folder).


WHAT IS SUMMARIZED IN THE SDD FILE SUMMARY:
All important header fields are interpreted and summarized in a text format at the beginning of the summary file.
The header summary is separated into 3 main subsections:
1. Incident Radiation Information - describes the type of radiation particle simulated.
2. Radiation Target Information - describes the cell/nucleus model that was irradiated in the simulation.
3. DNA Damage Information - describes the damage definition for the appearance of Double Strand Breaks, and other exposure information.

The data block is summarized under the subsection 'Chromosome Damages', where the following information is summarized:
1. The number of cell exposures contained in the SDD file (number of cells irradiated).
2. The number of damage entries (data rows) per cell exposure.
3. The associated damages per chromosome per cell exposure (i.e. number base damages, single-strand breaks, and double-strand breaks)
4. Total number of base damages, single-strand breaks, and double-strand breaks over all chromosomes per cell exposure.
5. The total number of base damages, single-strand breaks, and double-strand breaks over all chromosomes over all cell exposures.

WHAT IS PLOTTED IN THE SDD FILE KARYOGRAM:
<p align = "center">
<img src="karyogram_damages_vs_dose.gif" alt="Secondary electrons from 6 MeV photons at different doses"><br>
<em> Increasing DNA damages as a function of dose resulting from secondary electrons of a primary 6 MeV photon beam.</em>
</p>

1. All chromosomes passed by the user are plotted with distinct colors for visual clarity.
2. Centromeric positions using the '--karyogram human' option are realistic and adjusted based on the chromosome being drawn, whereas 
non-human genomes will be drawn using the generic centromere locations at approximately 50% of the chromosomes length. 
3. Double-strand breaks are denoted using black lines spanning slightly more than the width of the chromosome to disitnguish from single-strand breaks.
4. Single-strand breaks are denoted using white lines on the chromosomes and the line lengths are scaled based on the number of single-strand breaks 
within a given damage site (i.e. a given SDD data row) (usually between 0-5 SSBs). 
5. Summary at the top to describe the dose or fluence in the SDD file, the cell cycle phase, and the number of single- and double-strand breaks.
6. Legend at the bottom to describe what each symbol signifies in the karyogram.

WHAT IS SUMMARIZED IN THE SDR FILE SUMMARY:
1. SDR version, author, associated SDD file that produced the SDR file from the MEDRAS-MC output.
2. Number of chromosomes listed and their respective sizes in mega base pairs.
3. Per Cell summary of the cell subheader (number of double-strand breaks and misrepairs).
4. Number of mutations per each type.

WHAT IS PLOTTED IN THE SDR FILE KARYOGRAM 
1. Summary at the top, including which cell the Karyogram is plotted for, and the number of each mutation type.
2. Chromosomes plotted with their associated homologs if specified, and structural rearrangements from the associated mutations. 
3. The colors are used to depict where the chromosomal rearrangements came from, the white sections indicate where the deletions occurred on the 
original chromosomes. Extrachromosomal DNA (ecDNA) is drawn as circles with their diameters scaled to the length of the deleted segment for visual
clarity.
4. A legend of symbols at the bottom of the karyogram to denote what all the symbols mean, particulary for inversions and ecDNA.
5. So far, drawing is supported for long deletions, balanced inversions, balanced translocations, ecDNA, deletion-inversions, deletion-translocations,
and deletion-insertions, chromoplexy and chromothripsis. 

<p align = "center">
<img src="./examples/exampleSDR_cell0_karyogram_human.png" alt="Different mutation representations from irradiated DNA"><br>
<em> Depicting common DNA structural variations that arise from ionizing radiation.</em>
</p>


FURTHER KARYOGRAM PLOTTING WORK:
1. Add a zoom in and out and a panning option to conserve image quality. 
2. Add SDR mutation detection and drawing capabilities for BFB cycle.
3. Add functionality to be able to receive two X chromosomes instead of a Y and X chromosome and be able to plot them on the karyogram. 
5. Convert this command-line program into a linked c++ library to be called by any user who needs its functionality.
... and much more.

APPENDIX:
Example SDR file containing mutation entries for each type, the exact same as exampleSDR.txt, minus all the comments:
**** NOTE: In reality, the order of the SDR data entries do not matter, but they are logically ordered in the exampleSDR.txt file for clarity****

SDR Version, SDRv2.1;
Author, Samy Boutros;
Associated SDD File, None - Example for LongDels BalInvs BalTrans ecDNA DelInvs MultiDels MultiecDNA DelTrans DelIns Chromoplexy Chromothripsis;
Intact Chromosome Sizes, 46, 263.2, 258.1, 212.4, 203.4, 192.2, 181.8, 168.8, 155.5, 149.4, 144.0, 142.9, 140.8, 121.3, 113.0, 106.9, 94.7, 83.5, 81.0, 68.0, 66.6, 50.0, 52.6, 263.2, 258.1, 212.4, 203.4, 192.2, 181.8, 168.8, 155.5, 149.4, 144.0, 142.9, 140.8, 121.3, 113.0, 106.9, 94.7, 83.5, 81.0, 68.0, 66.6, 50.0, 52.6, 61.6, 165.2,
***end of master header***

***subheader - cell0***
Cell ID, 0;
Mutated Chromosome Sizes, {2: 212.4, 3: 203.4, 5: 181.8, 6: 168.8, 7: 155.5, 8: 149.4, 9: 144.0, 10: 142.9, 11: 140.8, 12: 121.3, 13: 113.0, 14: 106.9, 15: 94.7, 16: 83.5, 17: 81.0, 18: 68.0, 19: 66.6, 20: 50.0, 21: 52.6, 22: 263.2, 23: 258.1, 24: 212.4, 25: 203.4, 26: 192.2, 27: 181.8, 28: 168.8, 29: 155.5, 30: 149.4, 31: 144.0, 32: 142.9, 33: 140.8, 34: 121.3, 35: 113.0, 36: 106.9, 37: 94.7, 38: 83.5, 39: 81.0, 40: 68.0, 41: 66.6, 42: 50.0, 43: 52.6, 44: 61.6, 45: 165.2, 46: 148.461676, 47: 43.738324, 48: 263.2, 49: 358.985622, 50: 162.314378};
Intact Strands ID, 43, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45;
Total DSB Count, ;
Total Misrepair Count, ;
Medras-MC Log, ;

***data - cell0***
0; 1; 1/0/263.2/1; 1;
0; 2; 2/0/258.1/1; 1;
0; 3; 3/0/212.4/1; 1;
0; 4; 4/0/203.4/1; 1;
0; 5; 5/0/192.2/1; 1;
0; 6; 6/0/181.8/1; 1;
0; 7; 7/0/168.8/1; 1;
0; 8; 8/0/155.5/1; 1;
0; 9; 9/0/149.4/1; 1;
0; 10; 10/0/144.0/1; 1;
0; 11; 11/0/142.9/1; 1;
0; 12; 12/0/140.8/1; 1;
0; 13; 13/0/121.3/1; 1;
0; 14; 14/0/113.0/1; 1;
0; 15; 15/0/106.9/1; 1;
0; 16; 16/0/94.7/1; 1;
0; 17; 17/0/83.5/1; 1;
0; 18; 18/0/81.0/1; 1;
0; 19; 19/0/68.0/1; 1;
0; 20; 20/0/66.6/1; 1;
0; 21; 21/0/50.0/1; 1;
0; 22; 22/0/52.6/1; 1;
0; 23; 23/0/263.2/1; 1;
0; 24; 24/0/258.1/1; 1;
0; 25; 25/0/212.4/1; 1;
0; 26; 26/0/203.4/1; 1;
0; 27; 27/0/192.2/1; 1;
0; 28; 28/0/181.8/1; 1;
0; 29; 29/0/168.8/1; 1;
0; 30; 30/0/155.5/1; 1;
0; 31; 31/0/149.4/1; 1;
0; 32; 32/0/144.0/1; 1;
0; 33; 33/0/142.9/1; 1;
0; 34; 34/0/140.8/1; 1;
0; 35; 35/0/121.3/1; 1;
0; 36; 36/0/113.0/1; 1;
0; 37; 37/0/106.9/1; 1;
0; 38; 38/0/94.7/1; 1;
0; 39; 39/0/83.5/1; 1;
0; 40; 40/0/81.0/1; 1;
0; 41; 41/0/68.0/1; 1;
0; 42; 42/0/66.6/1; 1;
0; 43; 43/0/50.0/1; 1;
0; 44; 44/0/52.6/1; 1;
0; 45; 45/0/61.6/1; 1;
0; 46; 46/0/165.2/1; 1;

------------------ Multiple deletion per chromosome ---------------------
0; 47; 5/0/55.0/1, 5/88.0/137.4/0, 5/181.2/192.2/0; 1;
0; 48; 5/137.4/181.2/0; 1;
0; 49; 5/55.0/88.0/0; 1;
-------------------------------------------------------------------------

------------ Single and multiple inversions per chromosome --------------
0; 50; 1/0/90.203953/0, 1/110.236709/90.203953/0, 1/110.236709/263.2/1; 1;

0; 51; 3/0/100.8/1, 3/120.5/100.8/0, 3/120.5/150.9/0, 3/173.2/150.9/0, 3/173.2/212.4/0; 1;
-------------------------------------------------------------------------

----------------- Balanced translocation examples -----------------------
-- Bottoms Swapped:
0; 52; 2/0/223.579143/1, 1/127.793521/263.2/0; 1;
0; 53; 1/0/127.793521/1, 2/223.579143/258.1/0; 1;

-- Tops swapped:
0; 66; 30/0/60.0/1, 26/75.3/203.4/0; 1;
0; 67; 26/0/75.3/1, 30/60.0/155.5/0; 1;

-- Bottom swaps with top:
0; 68; 13/121.3/73.2/0, 9/20.7/149.4/1; 1;
0; 69; 13/0/73.2/1, 9/0/20.7/0; 1;

-- Central fragments swap:
0; 70; 12/0/70.0/1, 16/60.0/78.0/0, 12/100.0/140.8/0; 1;
0; 71; 16/0/60.0/1, 12/70.0/100.0/0, 16/78.0/94.7/0; 1;
-------------------------------------------------------------------------

---------------------------- ecDNA examples -----------------------------
-- Multiple deletions forming multiple ecDNA fragments:
0; 54; 36/0/55.0/1, 36/70.0/88.0/0, 36/100.0/113.0/0; 1;
0; 55; 36/55.0/70.0/0; 0;
0; 56; 36/88.0/100.0/0; 0;

-- Multiple deletions recombining to form one ecDNA fragment:
0; 57; 15/0/30.0/1, 15/50.0/70.0/0, 15/90.0/106.9/0; 1;
0; 58; 15/30.0/50.0/0, 15/70.0/90.0/0; 0;
-------------------------------------------------------------------------

-------------------- Deletion-Inversion ---------------------------------
0; 59; 8/0/44.0/0, 8/112.5/85.0/1, 8/112.5/155.5/0; 1;
0; 60; 8/44.0/85.0/1; 1;

----------------------- Deletion-Translocation --------------------------
0; 61; 10/0/65.0/1, 11/79.0/142.9/0; 1;
0; 62; 11/0/79.0/1, 10/100.0/144.0/0; 1;
0; 63; 10/65.0/100.0/0; 1;

----------------------- Deletion-Insertion ------------------------------
0; 64; 24/0/121.265341/1, 24/188.892635/258.1/0; 1;
0; 65; 23/0/130.124698/1, 24/121.265341/188.892635/0, 23/130.124698/263.2/0; 1;

-------------------------- Chromoplexy ----------------------------------
0; 72; 41/0/20.0/0, 39/15/83.5/1; 1;
0; 73; 40/0/70.0/1, 39/15.0/0/0; 1;
0; 74; 40/81.0/70.0/0, 41/20.0/68.0/1; 1;

-------------------------- Chromothripsis -------------------------------
-- Single-chromosome chromothripsis:
0; 75; 6/0/15/0, 6/25/15/0, 6/30/35/0, 6/50/62/1, 6/70/62/0, 6/70/75/0, 6/85/80/0, 6/98/105/0, 6/126/105/0, 6/135/140/0, 6/155/168/0, 6/174/181.8/0; 1;
0; 76; 6/25/30/0; 1;
0; 77; 6/35/50/0; 0;
0; 78; 6/75/80/0; 1;
0; 79; 6/85/98/0; 1;
0; 80; 6/126/135/0; 1;
0; 81; 6/140/155/0; 1;
0; 82; 6/168/174/0; 1;

-- 2-chromosome chromothripsis
0; 83; 25/0/35/0, 29/145/155/0, 25/50/65/0, 25/80/65/0, 25/85/95/1, 29/32/16/0, 25/130/200/0, 29/16/0/0; 1;
0; 84; 25/212.4/200/0, 29/32/60/1, 25/35/50/0, 29/130/70/0, 29/130/145/0, 29/155/168.8/0; 1;
0; 85; 25/80/85/0; 1;
0; 86; 25/98/130/0; 1;
0; 87; 29/60/70/0; 1;
-------------------------------------------------------------------------
