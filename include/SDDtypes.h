#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

// The SDD header may contain any of the following header fields.
struct Header
{
    std::string sdd_version;
    std::string software;
    std::string author;
    std::string simulation_details; 
    std::string source;	
    int source_type; 
    std::vector<int> incident_particles; 
    std::vector<double> mean_particle_energy; 
    std::string energy_distribution; 
    std::vector<double> particle_fraction; 
    std::vector<double> dose_or_fluence; 
    double dose_rate = 0.0;
    std::string irradiation_target;
    std::vector<double> volumes; 
    std::vector<double> chromosome_sizes;
    double DNA_density; 
    std::vector<double> cell_cycle_phase; 
    std::vector<int> DNA_structure; 
    int in_vitro_or_in_vivo; 
    int proliferation_status;
    std::vector<double> microenvironment;
    std::vector<double> damage_definition;
    double time;
    std::vector<int> damage_and_primary_count;
    std::vector<int> data_entries;
    std::string additional_information;

};

//--------------------------------------------------
// Data Entry 1
// Classification
//--------------------------------------------------
struct Classification
{
    int exposureMarker = 0;
    int eventID = 0;
};

//--------------------------------------------------
// Data Entry 2
// Spatial Position
//--------------------------------------------------
struct Position
{
    double x = 0.0;
    double x_min = 0.0;
    double x_max = 0.0;
    double y = 0.0;
    double y_min = 0.0;
    double y_max = 0.0;
    double z = 0.0;
    double z_min = 0.0;
    double z_max = 0.0;
};

//--------------------------------------------------
// Data Entry 3
// Chromosome IDs
//--------------------------------------------------
struct ChromosomeID
{
    int dnaStructure = 0;
    int chromosomeNumber = 0;
    int chromatidNumber = 0;
    int chromosomeArm = 0;
};

//--------------------------------------------------
// Data Entry 4
// Chromosome Position
//--------------------------------------------------
struct ChromosomePosition
{
    double position = 0.0; 
    bool isFractional = false; 
};

//--------------------------------------------------
// Data Entry 5
// Cause
//--------------------------------------------------
struct DamageCause
{
    int cause = 0;
    int numDirectDamages = 0;
    int numIndirectDamages = 0;
};

//--------------------------------------------------
// Data Entry 6
// Damage Types
//--------------------------------------------------
struct DamageType
{
    int numBaseDamages = 0;
    int numSingleBackboneBreaks = 0;
    int presenceOfDSB = 0;
};

//--------------------------------------------------
// Data Entry 7
// Full break spec
//--------------------------------------------------
struct FullBreakSpec
{
    // 1 = 5' to 3' backbone
    // 2 = 5' to 3' bases
    // 3 = 3' to 5' bases
    // 4 = 3' to 5' backbone
    std::vector<int> strand;

    // Identifies the DNA base where damage occurs
    std::vector<int> base;

    // 0 = no damage
    // 1 = direct damage
    // 2 = indirect damage
    // 3 = direct + indirect damage
    std::vector<int> baseDamageType;
};


// Exposure entry structure, need to include optional fields later.
struct DamageEntry
{
    Classification classification;
    Position position;
    ChromosomeID chromosomeID;
    ChromosomePosition chromosomePosition;
    DamageCause damageCause;
    DamageType damageType;
    FullBreakSpec fullBreakSpec;

    std::string rawLine;

};

struct Exposure
{
    int exposureID = 0;
    std::vector<DamageEntry> damages;
};


//--------------------------------------------------
// Chromosome damage summary
//--------------------------------------------------
struct ChromosomeDamageSummary
{
    int dnaStructure = 0;
    int chromosomeNumber = 0;

    int numBaseDamages = 0;
    int numSingleStrandBreaks = 0;
    int numDSBs = 0;
};
