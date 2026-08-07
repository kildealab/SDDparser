#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

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
    int exposureID = 0;
    int eventID = 0;
};

//--------------------------------------------------
// Data Entry 2
// Spatial Position
//--------------------------------------------------
struct Position
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

//--------------------------------------------------
// Data Entry 3
// Chromosome Information
//--------------------------------------------------
struct Chromosome
{
    int chromosomeID = 0;
    int chromatidID = 0;
    int chromosomeNumber = 0;
    int chromatidNumber = 0;
};

//--------------------------------------------------
// Data Entry 4
// Damage Cause
//--------------------------------------------------
struct DamageCause
{
    uint64_t cause = 0;
};

//--------------------------------------------------
// Data Entry 5
// DNA Segment
//--------------------------------------------------
struct DNASegment
{
    int chromosome = 0;
    int copy = 0;
    int arm = 0;
};

//--------------------------------------------------
// Data Entry 6
// Chromatid Segment
//--------------------------------------------------
struct ChromatidSegment
{
    int chromatid = 0;
    int copy = 0;
    int arm = 0;
};

//--------------------------------------------------
// Data Entry 7
// Damage Description
//--------------------------------------------------
struct DamageType
{
    int baseDamage = 0;
    int singleStrandBreak = 0;
    int doubleStrandBreak = 0;
};

struct DamageEntry
{
    Classification classification;
    Position position;
    Chromosome chromosome;
    DamageCause cause;
    DNASegment dnaSegment;
    ChromatidSegment chromatidSegment;
    DamageType damage;

    std::string rawLine;

};

struct Exposure
{
    int exposureID = 0;
    std::vector<DamageEntry> damages;
};
