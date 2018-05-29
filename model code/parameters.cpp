#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <random>
#include "Vector.hpp"

//#include <boost/math/constants/constants.hpp>
#include <cmath>

// temp run parms
//#include "atdparms.hpp"


//TODO: This typedef is unnecessary. The code in main.cpp should be updated to
//use Vector straight. 
typedef Vector vec_T;
const float_T pi= 3.14159265358979323846;
//const float_T pi = M_PI;


//This will define our centrosome options. 
enum MTOC {M_CENTROSOME, D_CENTROSOME};

//Debug parameters. 
bool spitValues = false;

//Type Parameters
const bool motherSpringOn = false;
const bool daughterSpringOn = false;
const bool translation = false;
const bool ONLY_COMMA = true;

//Parameters
const float_T Duration = 800;         //duration in minutes
const float_T Tau      = 1.0/4000.0; //time step in minutes

// MT Parameters
#define MT_numb_M 1157 //for wider upper env
//#define MT_numb_M 1000
#define MT_numb_D 1000

vec_T MT_Pos_M[MT_numb_M];
vec_T MT_Pos_D[MT_numb_D];
bool MT_Growing_M[MT_numb_M];
bool MT_Growing_D[MT_numb_D];
float_T MT_GrowthVel_M[MT_numb_M];
float_T MT_GrowthVel_D[MT_numb_D];



//  Contact Parameters:
//float_T contact_length = 400*Tau; //This quantity is measured in min
float_T contact_length = 0.1450;
//float_T contact_length = Duration;

//float_T contact_length_dynein = 400*Tau;// This quantity measures the contact time for each MT pullig due to dynein 1/koff
float_T contact_length_dynein = 0.06784;
//float_T contact_length_dynein = Duration;

//float_T contact_length_push = 30*Tau; // This quantity measures the contact time for each MT pushing on the cortex =1/fcat.
float_T contact_length_push = 0.00776;
//float_T contact_length_push = 1.5*Tau;

float_T MT_Contact_M[MT_numb_M];
float_T MT_Contact_D[MT_numb_D];

float_T MT_Contact_M_push[MT_numb_M];
float_T MT_Contact_D_push[MT_numb_D];

//  MT Growth/Shrinking Parameters from the DifferentSpringsAndRateTables.pdf
//const float_T Vg   = 0.5*60; //growth velocity in mum/min
const float_T Vg   = 44.08139;
//const float_T Vs_c = 2*60; //shortening velocity after contact in mum/min
const float_T Vs_c = 155.80173;
//const float_T Vs   = 0.25*60; //shortening velocity in mum/min
const float_T Vs   = 9.81674;
//const float_T kc   = .04*60;   //catastrophe frequency in min-1
const float_T kc   = 2.12004;
//const float_T kr   = .1*60;   //rescue frequency in min-1
const float_T kr   =  3.77174;

float_T Pr_catastrophe = 1-exp(-kc*Tau);
float_T Pr_rescue      = 1-exp(-kr*Tau);


//   Force Parameters
const float_T F_MT_pull   = 1.0; //Force per MT in pN.
//const float_T F_MT_pull   = 0.82545;
const float_T F_MT_push   = 0; //Force per MT in pN.
//const float_T F_MT_push   = 0.44531;

const float_T Fratio = 1.0;// Env. killing: 0.32299363312512985;
vec_T force_M;
vec_T force_D;
vec_T force;
float_T torque_M;
float_T torque_D;
float_T torque;


//   Envelope Parameters
//const float_T envWidthM = 2*pi/3;
const float_T envWidthM = 2.32404; //for wider upper env
const float_T envelopeM[2] = {pi/2.0 - envWidthM/2.0, pi/2.0 + envWidthM/2.0}; // The envelope in which MTs from M can grow. 
const float_T envWidthD = 2*pi/3;
const float_T envelopeD[2] = {pi/2.0 - envWidthD/2.0, pi/2.0 + envWidthD/2.0}; // The envelope in which MTs from M can grow.

//   Spring Parameters
const float_T kM = 3;
const float_T kD = 3;

// Pronucleus Parameters
const float_T R1_max = 16;       //Embryo width in mum
const float_T R2_max = 15;       //Embryo length in mum
const float_T Prad   = 4;        //Pronucleous radius (mum)
//const float_T Eta    = 1.25;     //translational drag Cytoplasmic viscosity =1/60 (pN min/min), Drag coeff = 6 pi Prad/60 when springs are off
const float_T Eta    = 6.60540;
//const float_T Mu     = 5*Eta;      //Rotational drag coeff ((pN mum)/min)
//const float_T Mu     = 10000;
const float_T Mu     = 532.28096;
const float_T Eta2   = Eta/10;   //Trans. drag coeff of each MTOC (pN min/mum) WHEN SPRINGS ARE ON
//const float_T Eta2   = __ATD_PARM_2__;
//const float_T Eta2   = 100;
const float_T kbT    = .00414;   //pN mum
const float_T D      = kbT/Eta2; //Diffusion Coefficient for each MTOC free motion. (mum^2/min)
const float_T Dp     = 100*kbT/Eta; //Diffusion Coefficient for the pronucleus
const float_T Dpr    = 100*kbT/Mu; //Roational Diffusion coefficient for the pronucleus

//  Starting Centered Coordinates:
//const float_T startPsi = pi/2.0;
//const float_T startX   = 0;
//const float_T startY   = 0;

//   Off-center Coordinates
const float_T startPsi = pi/2.0;
const float_T startX   = 0.4*R1_max; //at 70:30 mark
const float_T startY   = 0;

//  General Coordinate Initializations: 
float_T psi;
vec_T basePosM;
vec_T basePosD;
vec_T proNucPos;

// Band parameters: 
//  regionAngles:
//   The regionAngles vector defines the endpoints of the angular span of
//   various regions along the cortex. It MUST start with 0 and it MUST end with
//   2*pi. Each region is presumed to be homogenous and all encompassing: In
//   particular, a single band of a protein on the cortex does not totally
//   define a region unless it is the only protein defining connectivity
//   probability in that region. If two or more bands intersect, you must make
//   your regions various homogenous intersections such that the probability is
//   constant throughout any given region. Note that technically speaking, in
//   implementation a region is defined by two endpoints \theta_1, \theta_2 \in
//   [0,2\pi] such that a point on the cortex at angle \alpha is in that region
//   if \theta_1 <= \alpha < \theta_2
//  regionProbabilities: 
//    This defines the probabilities associated with the regions defined by the
//    regionAngles variable. It has length one less than the regionAngles
//    vector, as it is broken up into regions, not enpoints of regions. 
//  regionForceMultipliers:
//    This defines the force multipliers associated with the regions defined by
//    the regionAngles variable. It has length one less than the regionAngles
//    vector, as it is broken up into regions, not enpoints of regions. 

//No Bands:
//const int numRegions = 1;
//const float_T regionAngles[numRegions + 1] = {0, 2*pi};
//const float_T regionProbabilities[numRegions] = {1};
//const float_T regionForceMultipliers[numRegions] = {1};

//Anterior to posterior difference only (no push bands):
const int numRegions = 3;
//R1_max = 16, R2_max = 15;
const float_T myPos = 1.3564; //this is the 60:40 cartesian angle
const float_T regionAngles[numRegions+1] = {0, myPos, 2*pi - myPos, 2*pi};
const float_T regionProbabilities[numRegions] = {1,0.65,1}; //lower ant binding probability
const float_T regionForceMultipliers[numRegions] = {1,1,1}; //all pulling

//Standard Bands:
//const int numRegions = 5;
////R1_max = 16, R2_max = 15;
//const float_T start = 1.2439; //x=5.138
//const float_T end = 1.4946; //x=1.218
//const float_T regionAngles[numRegions+1] = {0, start, end, 2*pi - end, 2*pi - start, 2*pi}; //push bands
//const float_T regionProbabilities[numRegions] = {1,1,0.65,1,1}; //lower ant binding probability
//const float_T regionProbabilities[numRegions] = {1,1,1,1,1}; //no diff in binding probability
//const float_T regionForceMultipliers[numRegions] = {1,-1,1,-1,1}; //for push band force multipliers

// MT density limitations
//  Only one contact per window, windows of length ~1
//R1_max = 16, R2_max = 15;
const size_t numberContactWindows = 96;
const float_T contactWindowAngles[numberContactWindows+1] =
{0,0.0635,0.127,0.1905,0.254,0.318,0.382,0.446,0.5105,0.575,0.64,0.705,0.7705,0.836,0.902,0.968,1.0345,1.101,1.168,1.235,1.3025,1.37,1.4375,1.505,1.5725,1.64,1.7075,1.775,1.8425,1.91,1.977,2.044,2.1105,2.177,2.243,2.309,2.3745,2.44,2.505,2.57,2.6345,2.699,2.763,2.827,2.891,2.9545,3.018,3.0815,3.145,3.2085,3.272,3.3355,3.399,3.463,3.527,3.591,3.6555,3.72,3.785,3.85,3.9155,3.981,4.047,4.113,4.1795,4.246,4.313,4.38,4.4475,4.515,4.5825,4.65,4.7175,4.785,4.8525,4.92,4.9875,5.055,5.122,5.189,5.2555,5.322,5.388,5.454,5.5195,5.585,5.65,5.715,5.7795,5.844,5.908,5.972,6.036,6.0995,6.163,6.2265,6.283185};
bool contacts[numberContactWindows];

// File Parameters
std::ofstream file;
std::string fileName  = "";
std::string fileDir   = "../data/";
std::string fileOrder = "t,proNucPos,psi,MT_Pos_M,MT_Pos_D,force_M,force_D,\
                         force,torque_M,torque_D,torque,basePosM,basePosD";

// Random Number Generation Parameters: 
std::normal_distribution<float_T> stdNormalDist;
std::uniform_real_distribution<float_T> stdUniformDist(0.0,1.0);
std::default_random_engine generator;
