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
//#define MT_numb_M 1284
#define MT_numb_M 1157 //for wider upper env
//#define MT_numb_M 1000
//#define MT_numb_M 684
#define MT_numb_D 1000
//#define MT_numb_D 527

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
const float_T R1_max = 15;       //Embryo width in mum
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
const float_T startPsi = pi/2.0;
const float_T startX   = 0;
const float_T startY   = 0;

//   Off-center Coordinates
//const float_T startPsi = pi/2.0;
//const float_T startX   = 0.4*R1_max; //at 70:30 mark
//const float_T startY   = 0;

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
//R1_max = R2_max = 15
const float_T myPos = 1.3694; //this is the 60:40 cartesian angle
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
//R1_max = R2_max = 15
const size_t numberContactWindows = 94;
const float_T contactWindowAngles[numberContactWindows+1] =
{0,0.066842,0.13368,0.20053,0.26737,0.33421,0.40105,0.4679,0.53474,0.60158,0.66842,0.73527,0.80211,0.86895,0.93579,1.0026,1.0695,1.1363,1.2032,1.27,1.3368,1.4037,1.4705,1.5374,1.6042,1.6711,1.7379,1.8047,1.8716,1.9384,2.0053,2.0721,2.139,2.2058,2.2726,2.3395,2.4063,2.4732,2.54,2.6069,2.6737,2.7405,2.8074,2.8742,2.9411,3.0079,3.0748,3.1416,3.2084,3.2753,3.3421,3.409,3.4758,3.5426,3.6095,3.6763,3.7432,3.81,3.8769,3.9437,4.0105,4.0774,4.1442,4.2111,4.2779,4.3448,4.4116,4.4784,4.5453,4.6121,4.679,4.7458,4.8127,4.8795,4.9463,5.0132,5.08,5.1469,5.2137,5.2805,5.3474,5.4142,5.4811,5.5479,5.6148,5.6816,5.7484,5.8153,5.8821,5.949,6.0158,6.0827,6.1495,6.2163,6.283185};
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
