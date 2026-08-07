#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

struct Header
{
    std::string version;
    std::string software;
    std::string author;

    std::vector<double> chromosomeSizes;
    std::vector<int> dataEntries;

    double dose = 0.0;
    double doseRate = 0.0;

    std::string sourceType;
    std::string irradiationTarget;

    // ...add the remaining SDD header fields here...
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
