// Copyright Matus Krajnak 2018
// Pixelated DPC is distributed under the GNU General Public License (GPL)
// If you use Pixelated DPC, we kindly ask that you cite the following:

// 1. Pixelated detectors and improved efficiency for magnetic
// imaging in STEM differential phase contrast, M Krajnak, D McGrouther,
// D Maneuski, V O'Shea, S McVitie, Ultramicroscopy 165, 42-50 (2016)

// 2. Advanced detection in Lorentz microscopy: pixelated detection 
// in differential phase contrast scanning transmission electron 
// microscopy, M Krajnak, PhD Thesis, University of Glasgow (2017)


//main file calling important functions and calculation

#include "dpc_support.h"

using namespace std;
using namespace af;

//main function
int main(int argc, char** argv ) {
    //program info
    
    af::Window lala;
    showCreator();
    //create and check a struct s - everything can be passed to functions (defined in dpc_support.h/.cpp)
    DefaultVals s; 
    //set the important variables by passing onto function from dpc_support.cpp and check the basics for running the analysis
    setDefaultSettings(s);    
    //read setting from command line
    getSettings(s,argc,argv);

    // if no options given show help and exit
    if(argc == 1) {
        showhelpinfo();
        exit(1);
    }

    //check if it is good to go
    check_s(s);
    
    //tell what compute engine is selected by arrayfire
    cout << "\n"; af::info(); cout << "\n";

    
    //OPEN BIG BINARY FILE FOR ANALYSIS - open at the end for filelength check
    ifstream fs;
    fs.open(s.binfile.c_str(), ios::out | ios::binary | ios::ate);
    //file checks and go back to beginning
    data_check(fs, s);
        
    if (s.AUTO) {
        //generate mask file automatically
        s.mask = autoMask(fs, s);
    } else {
        //or get mask from a manually created file (using imageJ or sorts)
        //get mask file from maskfile - FreeImage supported formats
        s.mask = get_mask(s);
    }
    
    //prepare Edge = fft(edge) of mask for efficiency in cross-correlation
    s.Edge = mask_Edge(s);

    //create window for visualisation
    Window myWindow(256*4, 64 + 256*2, "computation visualisation");
    if (!s.show) myWindow.setVisibility(false);
    
    //initiate Results
    Results calcs;
    //generate check variable for exceptions
    bool done = false;

    //show that calculation has started
    calcString();
    
    //time operational part of the code
    s.start = timer::start();

    //while / try / catch is done to correct for running out of memory on gpu - it may not work 100%
    while (!done) {
        try {
            //go to the beginning of the file and restart frame counter
            fs.seekg(0, ios::beg);
            s.frame_i = 0;
            //arrays to store resulting calculations
            calcs = createResults(s);
            
            // ======CALCULATION START=======
            //       bulk calculation
            compute(fs, s, calcs, myWindow);
            //deal with any leftovers at the end of the file if XY % STEP > 0
            if (s.leftovers > 0) {
                s.STEP = s.leftovers;
                compute(fs, s, calcs, myWindow, true);
            }
            // ======CALCULATION FINISH=======
            done = true;
            
        // if there is no memory on GPU exception this will halve s.STEP and try again
        } catch (af::exception& ae) {
            if ((ae.err() == 101 || 998) and (s.STEP > 1)) {
                s.STEP=int(s.STEP/2);
                s.leftovers = s.SCAN_XY % s.STEP;
                af::deviceGC();
                cout << "\nAF error 101 or 998: device run out of memory\n";
                cout << "Trying with parralel frames halved to: " << s.STEP << endl << endl;
            } else {
                //other error should just complain and not run compute
                std::cerr << ae.what() << std::endl;
                done = true;
            }
        }
    }

    float finito = timer::stop();
    cout << endl;
    printf("\relapsed seconds: %g \n", roundFloat(finito));
    cout << "average fps " << roundFloat(s.SCAN_XY / finito) << endl<<endl;

    //check if the directory exists and if yes add 0 at the end to avoid the loss of data. When created, chdir into it
    while (DirectoryExists(s.dirname.c_str())) s.dirname+="0";
	mkdir(s.dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    chdir(s.dirname.c_str());
    
    //save and show resulting images
    saveimages(s,calcs);
    af::deviceGC();
    af::sync();
    //close the binary file
    fs.close();
    
    return 0;
}
