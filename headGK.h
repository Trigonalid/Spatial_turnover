#ifndef MMHEADPCOL_H
#define MMHEADPCOL_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <numeric>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

const double PI = 3.14159265;
const double TWOPI = 6.28318531;

// Stores total colonization and extinction rates
struct totals {
    double Hcol; // Total host colonization rate
    double Pcol; // Total parasitoid colonization rate
    double Hext; // Total host extinction rate
    double Pext; // Total parasitoid extinction rate
};

// Represents a patch in the metacommunity
struct Patch {
    double x, y; // Spatial coordinates
    vector<int> Hpres; // Host presence (1=present, 0=absent)
    vector<int> Ppres; // Parasitoid presence
    vector<double> Hcol; // Host colonization rates
    vector<int> Pcol; // Parasitoid colonization counts
    vector<double> PcolDub; // Parasitoid colonization rates
    vector<double> Hext; // Host extinction rates
    vector<int> Pext; // Parasitoid extinction counts
    vector<double> PextDub; // Parasitoid extinction rates
    vector<vector<int>> HconnecIND; // Host connectivity indices
    vector<vector<int>> PconnecIND; // Parasitoid connectivity indices
};

// Simulation timing parameters
struct Times {
    int interval; // Time step interval
    double maxT; // Maximum simulation time
    int rec; // Recording interval
    double totalruntime; // Maximum runtime in seconds
    int N; // Number of replicates
};

// Model parameters
struct Pars {
    int H, P; // Number of host and parasitoid species
    vector<double> mH; // Host colonization rates
    vector<double> mP; // Parasitoid colonization rates
    vector<double> xH; // Host extinction rates
    vector<double> xP; // Parasitoid extinction rates
    vector<double> LP; // Parasitoid dispersal distances
    vector<double> LH; // Host dispersal distances
    double U; // Spatial domain size
    double gamma; // Effect of hosts on parasitoid extinction
    double delta; // Effect of parasitoids on host extinction
    int numvil; // Number of patches
    int set; // Simulation set identifier
    int pcol; // Parasitoid colonization mode (0=fixed, 1=density-dependent)
    vector<vector<int>> ParaHosts; // Hosts for each parasitoid
    vector<vector<int>> HostPara; // Parasitoids for each host
    string hostfile, parafile; // Input file names
};

// Function declarations
double RunOnceInt(double);
double random_exp(double);
void record(int);
double distance(double, double, double, double, double);
void RunMaxT(int);
void RunNReps(int);
void Initiate();
void PColonisation(int, int);
void HColonisation(int, int);
void PExtinguish(int, int);
void HExtinguish(int, int);
void UpdatePext();
void UpdatePcol();
double OneStep();
int IRandom(int, int);
double Random();
double chop(double);
int pickfromvec(const vector<double>&);

#endif
