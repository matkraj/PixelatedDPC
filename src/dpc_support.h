// Copyright Matus Krajnak 2018
// Pixelated DPC is distributed under the GNU General Public License (GPL)
// If you use Pixelated DPC, we kindly ask that you cite the following:

// 1. Pixelated detectors and improved efficiency for magnetic
// imaging in STEM differential phase contrast, M Krajnak, D McGrouther,
// D Maneuski, V O'Shea, S McVitie, Ultramicroscopy 165, 42-50 (2016)

// 2. Advanced detection in Lorentz microscopy: pixelated detection 
// in differential phase contrast scanning transmission electron 
// microscopy, M Krajnak, PhD Thesis, University of Glasgow (2017)


//definitions of functions

#ifndef DPC_SUPPORT_H
#define DPC_SUPPORT_H

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <arrayfire.h>
#include <random>
#include <dirent.h>
#include <cstring>

using namespace std;

//struct to store all analysis settings and arrays for comparison
struct DefaultVals{
    std::string dirname,binfile,maskfile,tmp11,tmp22,interactiveDOG;
    bool DOG,endian,sign,checkDOG1,checkDOG2,show,AUTO,hann;	 	
    unsigned int HEADERSIZE,FOOTERSIZE;	     
    int STEP,ind;		
    unsigned long ZEROPAD, SIZE_DET, SCAN_X, SCAN_Y, SCAN_XY, BYTESIZE_DET, frame_i,leftovers;	
    double sigma,dogSmall,dogLarge,DOWNSCALE,smooth_precise;
    uint bitdepth,bytedepth;
    af::timer start;
    af::array mask,edge,Edge;
    
    //AUTO settings
    unsigned int number_edges, edge_select, number_smooths, smooth_select, summing;
};

//struct to store resulting images
struct Results{
    af::array maxs,subs,sums,poss,dpcx,dpcy;
};

//some text showing functions
void showhelpinfo();
void showCreator();
void calcString();
void calcFinitoString();

//some helper functions
float roundFloat(float in);
std::string sConvert (float number);
bool DirectoryExists(const char* pzPath );
af::array smoothed_disk(int size, int radius, float smooth);
af::array zoom_image(const af::array &disk, float zoom);
af::array generate_hann(const af::array &disk);
af::array zeropad(const DefaultVals &s, const af::array &ar_line_small);
void endianSwap(const DefaultVals &s, int LENGTH, char* tmp);
af::array pointerToArray(DefaultVals s, char* tmpNoHead, bool mask = false);
void visualisation_window(const af::array &a, const af::array &b, const af::array &c, const af::array &d, const af::array &e, const af::array &f, const af::array &g, const af::array &h, af::Window& myWindow);
void saveimages(const DefaultVals &s, Results calcs);
af::array norm_hist(const af::array &in, float percentage);
void plot1d(af::Window &win, const af::array &ar1d, std::string lala = "1d data");
af::array downscale(const DefaultVals &s,const af::array &ar_line_small);
void gotHere();
void visualise_calc(DefaultVals &s, Results &results, af::Window& myWindow, const af::array &disk, const af::array &edge, const af::array &ccor, long indL, long indH);


// function acting on simulation setup
void data_check(std::ifstream& fs, const DefaultVals &s);
void generate_outputs(const DefaultVals &s,Results &calcs);
Results createResults(const DefaultVals &s);
void createLogFile(const DefaultVals &s);
void check_s(DefaultVals &s);
af::array get_mask(const DefaultVals &s);
void setDefaultSettings(DefaultVals &s);
int getSettings(DefaultVals &s, int argc, char** argv);

//image normalisations
af::array normalizeImage(const af::array &in);
void normalizeImageInPlace(af::array &in);
af::array normalizeCC(const af::array &in);
void normalizeCCinPlace(af::array &in);
void normaliseViewInPlace(af::array& disk);
af::array normaliseView(const af::array& disk);
af::array maxnorm(const af::array &disk);

//correlation functions
af::array partialConvolve(const DefaultVals &s, const af::array &edge_line);
af::array cross_correlate(const af::array &a1, const af::array &a2);
af::array cross_correlate2(const af::array &a1, const af::array &a2);
af::array smooth(const af::array &in, float smooth = 1.0, af::array smooths = af::constant(1,1));
af::array smooth1D(const af::array &in, float smooth = 1.0);
af::array generate_edge(DefaultVals s, af::array mask);
af::array generate_edge(const af::array &disk, float smooth = 1.0);
af::array mask_Edge(DefaultVals &s);

//reading binary data functions
af::array readDiskLine(std::ifstream& fs, const DefaultVals &s);
af::array get_disks_RANDOM(std::ifstream& fs,DefaultVals s, unsigned long noFrames);

//edge autodetection functions
af::array removeHotPixels(const af::array &in, float thr);
af::array select_disk(std::ifstream &fs, const DefaultVals &s, af::Window &myWindow, const af::array &disks, const af::array &edges);
af::array thr_range(const af::array& in, unsigned long bins);
void generate_disks_edges(std::ifstream &fs, const DefaultVals &s,af::Window &myWindow,af::array &disks,af::array &edges);
float find_smooth(std::ifstream &fs, const DefaultVals &s,af::Window &myWindow, const af::array &disk0, float refine = 0.0);
af::array finalise_disk(const DefaultVals &s, const af::array &disk0);
af::array autoMask(std::ifstream &fs, DefaultVals &s);

//main compute function
int compute(std::ifstream& fs, DefaultVals &s, Results &results, af::Window& myWindow, bool end = false);

    

#endif
