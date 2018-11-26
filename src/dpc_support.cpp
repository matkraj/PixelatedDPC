// Copyright Matus Krajnak 2018
// Pixelated DPC is distributed under the GNU General Public License (GPL)
// If you use Pixelated DPC, we kindly ask that you cite the following:

// 1. Pixelated detectors and improved efficiency for magnetic
// imaging in STEM differential phase contrast, M Krajnak, D McGrouther,
// D Maneuski, V O'Shea, S McVitie, Ultramicroscopy 165, 42-50 (2016)

// 2. Advanced detection in Lorentz microscopy: pixelated detection 
// in differential phase contrast scanning transmission electron 
// microscopy, M Krajnak, PhD Thesis, University of Glasgow (2017)


//this file contains support and computing function for pixDPC

#include "dpc_support.h"

using namespace std;
using namespace af;


// help output with a couple of examples
void showhelpinfo(){
  cout<<"Usage:   "<<" [-option] [argument] (default)"<<endl;
  cout<<"option:  "<<"-h  show help information"<<endl<<endl;
  cout<<"         "<<"-o  [DIR], output dir (outputDir)"<<endl;
  cout<<"         "<<"-a  [8/16/32], bitdepth of binary stream (16)"<<endl;
  cout<<"         "<<"-b  <FILE>, binary stream file"<<endl;
  cout<<"         "<<"-w  [0/1], endian swap - leave as is 0, swap 1, (1)"<<endl;
  cout<<"         "<<"-l  [0/1], unsigned/signed switch - unsigned 0, signed 1 (1)"<<endl;
  cout<<"         "<<"-A  auto generate mask - additional settings bellow"<<endl; 
  cout<<"         "<<"\t-S  [N], number of frames used for average (32)"<<endl;
  cout<<"         "<<"\t-N  [N], number of edges generated for statistical comparison (64)"<<endl;
  cout<<"         "<<"\t-M  [N], number of votes to select an edge (32)"<<endl;
  cout<<"         "<<"\t-O  [N], number of different smooth settings evenly spread to explore (32)"<<endl;
  cout<<"         "<<"\t-T  [N], number of votes to select a smooth value (32)"<<endl;
  cout<<"         "<<"\t-F  test with fast settings (-S 4, -N 4, -M 4, -O 4, -T 4)"<<endl;  
  cout<<"         "<<"-m  [FILE], mask file (will be ignored if -A)"<<endl;
  cout<<"         "<<"-e  [N], size of header of detector frame in bytes in front of each frame (384)"<<endl;
  cout<<"         "<<"-f  [N], size of footer of detector frame in bytes in the end of each frame (0)"<<endl<<endl; 
  cout<<"         "<<"========EDGE FILTER SETTINGS======="<<endl;
  cout<<"         "<<"-s  [f], sigma of Gaussian (0.0)"<<endl;
  cout<<"         "<<"-g  use DOG filter, -t and -u to specify smooth"<<endl;
  cout<<"         "<<"\t-t  [f], DOG sigma1 (0.5)"<<endl;
  cout<<"         "<<"\t-u  [f], DOG sigma2 (1.0)"<<endl<<endl;
  cout<<"         "<<"=======DETECTOR/STEM SCAN SETTINGS======"<<endl;
  cout<<"         "<<"-p  [N], number of frames processed in parallel (128)"<<endl;
  cout<<"         "<<"-x  [N], 'x' dimension of STEM scan (257)"<<endl;
  cout<<"         "<<"-y  [N], 'y' dimension of STEM scan (256)"<<endl;
  cout<<"         "<<"-d  [N], size of square detector (256)"<<endl;
  cout<<"         "<<"-z  [N], zero padding in multiples of detector size -d (2)"<<endl;
  cout<<"         "<<"-c  [f], downscale, use 2 to halve the detector size for speed -c >=1.0 (1.0)"<<endl<<endl;


  cout<<endl<<"         ===========Minimal example usage (Medipix3, 257x256 scan)===========\n         $ ./pixDPC -b bin_file"<<endl;
}

float roundFloat(float in) {
    long inn = std::round(in*100);
    return inn/100.0;
}

af::array normalizeImage(const af::array &in) {
    float min = af::min<float>(in);
    float max = af::max<float>(in);
    return 65635.0f*((in - min) / (max - min));
}

void normalizeImageInPlace(af::array &in) {
    float min = af::min<float>(in);
    float max = af::max<float>(in);
    in = 65635.0f*((in - min) / (max - min));
}

af::array normalizeCC(const af::array &in) {
    unsigned int x = in.dims(0);
    unsigned int y = in.dims(1);
    unsigned int z = in.dims(2);
    // make2D array
    af::array in2d = af::moddims(in,x*y,z);
    af::array mean2d = af::mean(in2d,0);
    af::array stdev2d = af::stdev(in2d,0);
    in2d = ( in2d - af::tile(mean2d,x*y) ) / af::tile(stdev2d,x*y);
    return af::moddims(in2d,x,y,z);
}

void normalizeCCinPlace(af::array &in) {
    unsigned int x = in.dims(0);
    unsigned int y = in.dims(1);
    unsigned int z = in.dims(2);
    in = af::moddims(in,x*y,z);
    af::array mean2d = af::mean(in,0);
    af::array stdev2d = af::stdev(in,0);
    in = ( in - af::tile(mean2d,x*y) ) / af::tile(stdev2d,x*y);
    in = af::moddims(in,x,y,z);
}

af::array generate_hann(const af::array &disk) {
    //this whole thing is probably wrong for 2D images
    af::array rangeX = af::range(disk.dims(0));
    af::array rangeY = af::range(disk.dims(1));

    rangeX = 0.5*(1-af::cos(2*M_PI*maxnorm(rangeX)));
    rangeY = 0.5*(1-af::cos(2*M_PI*maxnorm(rangeY)));
    
    rangeX = af::tile(rangeX,1,disk.dims(0));
    rangeY = af::reorder(rangeY,1,0);
    rangeY = af::tile(rangeY,disk.dims(1),1);
    
    af::array hann = rangeX * rangeY;
    hann = 1.0;
    return maxnorm(hann);
}

af::array norm_hist(const af::array &in, float percentage) {
    af::array out = in.copy();
    float mn = min<float>(out);
    float mx = max<float>(out);
    
    int noBins = 1024;
    af::array hist = histogram(out,noBins,mn,mx);
    
    int index = 0;
    int sumd = 0;
    float maxsum = out.dims(0) * out.dims(1) / 100 * percentage; // break for histogram stretch - to cover for errors in analysis
    
    //calculate minimum for resulting array
    while (sumd < maxsum) {
        sumd += sum<int>(hist(seq(index,index)));
        index++;
    }
    float minb = mn + index * ( (mx-mn) / noBins);

    index = 0;
    sumd = 0;
    //calculate maximum for resulting array
    while (sumd < maxsum) {
        sumd += sum<int>(hist(seq(noBins-index-1,noBins-index-1)));
        index++;
    }
    float maxb = mx - index * ( (mx-mn) / noBins);
    out = 65635.0f * ((out - minb) / (maxb - minb));
    return out;
}

af::array smoothed_disk(int size, int radius, float smooth) {
    if (radius == 0) {
        return af::gaussianKernel( size, size , smooth, smooth);
    } else {
        af::array x = range(size) - size/2;
        af::array y = af::reorder(x,1,0);
        
        x = af::tile(x,1,size);
        y = af::tile(y,size,1);
            
        af::array disk = (sqrt(x*x+y*y)<radius)*1.0f;
        af::array gaussian = af::gaussianKernel( size, size , smooth, smooth);
        disk = fftConvolve2(gaussian,disk,AF_CONV_DEFAULT);
        
        return maxnorm(disk);
    }
}

void normaliseViewInPlace(af::array& disk) {
    disk = disk - af::min<float>(disk);
    disk = disk / af::max<float>(disk);
}

af::array normaliseView(const af::array& disk) {
    af::array disk2 = disk - af::min<float>(disk);
    disk2 = disk2 / af::max<float>(disk2);
    return disk2;
}

af::array normaliseViewRange1D(const af::array& in, unsigned long maxrange) {
    int bottom = 0;
    if (maxrange>300) bottom = 260;
    af::array out = in - af::min<float>(in(seq(bottom,maxrange)));
    out = out / af::max<float>(out(seq(bottom,maxrange)));
    return out;
}

af::array smooth(const af::array &in, float smooth, af::array smooths) {
    af::array gaussian;
    if (smooths.dims(0) == 1) {
        gaussian = af::gaussianKernel(in.dims(0),in.dims(1), smooth, smooth);
        return af::fftConvolve2(in,gaussian,AF_CONV_DEFAULT);
    } else {
        af::array smooths3D = af::constant(0,in.dims(0),in.dims(1),smooths.dims(0));
        for (int x=0;x<smooths.dims(0);x++) {
            gaussian = af::gaussianKernel(
                in.dims(0),in.dims(1), 
                af::sum<float>(smooths(seq(x,x))), af::sum<float>(smooths(seq(x,x)))
            );
            smooths3D.slice(x) = af::fftConvolve2(in,gaussian,AF_CONV_DEFAULT);
        }
        return smooths3D;
    }
}
    
af::array smooth1D(const af::array &in, float smooth) {
    af::array gaussian = af::gaussianKernel(in.dims(0),in.dims(0), smooth, smooth);
    gaussian = gaussian(seq(),seq(gaussian.dims(0)/2,gaussian.dims(0)/2));
    return af::fftConvolve2(in,gaussian,AF_CONV_DEFAULT );
}

bool DirectoryExists(const char* pzPath) {
    if ( pzPath == NULL) return false;
    DIR *pDir;
    bool bExists = false;
    pDir = opendir (pzPath);
    if (pDir != NULL){
        bExists = true;
        (void) closedir (pDir);
    }
    return bExists;
}

af::array maxnorm(const af::array &disk) {
    return disk.as(f32)/af::max<float>(disk);
}

af::array max_thr_correlate(const af::array &in,const af::array &diskmed) {
    af::array maxs, disk_test, diskmed_edge, ccor;
    unsigned long bins = in.dims(2);
    unsigned long step = 32;
    unsigned long count = 0;
    
    while (count < bins) {
        //generate edges
        disk_test = generate_edge(in(seq(),seq(),seq(count,count+step-1)));
        count += step;
        diskmed_edge = generate_edge(diskmed);
        
        //normalise for cross-correlation disk - mean, disk / stdv
        normalizeCCinPlace(disk_test);
        normalizeCCinPlace(diskmed_edge);
        
        //tile edge generated from raw file
        diskmed_edge = af::tile(diskmed_edge,1,1,disk_test.dims(2));

        //cross-correlate
        ccor = cross_correlate2(disk_test,diskmed_edge);
        
        //scale by global maximum (just so the graph don't jump too much)

        //return maxima of cross-correlations - cross-corr spectrum (?)        
        maxs = af::join(0,maxs,af::flat(af::max(af::max(ccor,0),1)));
        af::deviceGC();
    }
    return maxnorm(maxs);
}

af::array thr_range(const af::array& in, unsigned long bins) {
    // create  <0,1> range
    af::array thresholds = af::range(bins)/(bins-1);

    // change to match the range of data
    float min_in = af::min<float>(in);
    float max_in = af::max<float>(in);
    thresholds = thresholds*(max_in - min_in) + min_in;
    // match the dimensions of the data
    thresholds = af::reorder(thresholds,1,2,0);
    thresholds = af::tile(thresholds,in.dims(0),in.dims(1));
    
    // copy the disk image over the range of thresholds for testing
    af::array disk_test = af::tile(in,1,1,thresholds.dims(2));
    
    // apply thresholds
    disk_test = 1.0f*(disk_test > thresholds);
    af::deviceGC();
    return disk_test;
}

unsigned long posMax1D(const af::array& in) {
    af::array maximum, position;
    af::max(maximum,position,af::flat(in));
    return af::sum<unsigned long>(position);
}

unsigned long posMin1D(const af::array& in) {
    af::array maximum, position;
    af::max(maximum,position,af::flat(in));
    return af::sum<unsigned long>(position);
}

af::array integral1D(const af::array& in) {
    af::array out = in.copy();
    for (int x=0;x<in.dims(0);x++)
        out(seq(x,x)) = af::sum<float>(in(seq(0,x)));
    return out;
}

void generate_disks_edges(ifstream &fs, const DefaultVals &s, af::Window &myWindow, af::array &disks, af::array &edges) {
    af::array disk,diskmed,hist,disk_test,diskmthr,maxs,integral_hist,maxs_filt,diskmthr_edge;
    int n = 0;
    while (!myWindow.close()) {
        //read disks
        disk = get_disks_RANDOM(fs,s,s.summing);
        //average disks
        disk = af::sum(disk,2);
        disk = removeHotPixels(disk,0.001);
        //filter for median - to remove any shot noise pixels
        diskmed = af::medfilt(disk);
        //show disk
        myWindow(0,0).image(normaliseView(diskmed),"disk median");
        //use hanning window
        //if (s.hann) disk = disk * generate_hann(disk);
        //select number of threshold bins
        unsigned int bins = 256;
        //compute histogram
        hist = af::histogram(diskmed,bins);
        // algorithm to create different threshold disks which according to equal distances in histogram
        disk_test = thr_range(diskmed,bins);
        // find the max of the cross-correlation for ideal threshold
        maxs = max_thr_correlate(disk_test,diskmed);
        //plot
        myWindow(0,1).hist(maxs,0,10,"thr_spectrum");
            
        //find position of maximum of the correlation spectrum
        unsigned long posmax = posMax1D(maxs);
        
        //apply the threshold to median filtered disk
        float threshold = (af::max<float>(diskmed) - af::min<float>(diskmed))*posmax/bins + af::min<float>(diskmed);
        diskmthr = 1.0f*(diskmed > threshold);
        //do another median for shot denoise at the edge
        diskmthr = af::medfilt(diskmthr,3,3);
        //apply edge filter
        diskmthr_edge = generate_edge(diskmthr,2.0);
        if (s.hann) diskmthr_edge = diskmthr_edge * generate_hann(diskmthr_edge);
        normalizeCCinPlace(diskmthr_edge);
        //add edge to edge images
        disks = af::join(2,disks,diskmthr);
        edges = af::join(2,edges,diskmthr_edge);

        normaliseViewInPlace(diskmthr);
        myWindow(0,2).image(diskmthr, "thr_disk");
        myWindow.show();
        
        af::sync();
        af::deviceGC();

        if (disks.dims(2)==s.number_edges) break;
    }
    af::deviceGC(); 
}

af::array select_disk(ifstream &fs, const DefaultVals &s,af::Window &myWindow, const af::array &disks, const af::array &edges) {
    af::array maxs = af::constant(0.0,edges.dims(2));
    int n = 0;
    af::array disk, d_edge, ccor, max_spectrum;
    while (!myWindow.close()) {
        disk = get_disks_RANDOM(fs,s,1);
        disk = af::medfilt(disk);
        d_edge = generate_edge(disk);
        d_edge = af::medfilt(d_edge);
        if (s.hann) d_edge = d_edge * generate_hann(d_edge);
        d_edge = normalizeCC(d_edge);
        myWindow(0,0).image(normaliseView(d_edge), "edge");
        
        d_edge = af::tile(d_edge,1,1,edges.dims(2));
        ccor = cross_correlate2(d_edge,edges);
        max_spectrum = af::flat(af::max(af::max(ccor,0),1));
        unsigned int mxpos = posMax1D(max_spectrum);
        
        maxs(seq(mxpos,mxpos))+=1;
        
        myWindow(0,1).hist(10*maxs,0,maxs.dims(0));
        myWindow(0,2).image(normaliseView(edges.slice(mxpos)),"best");
        myWindow.show();
        af::deviceGC();
        if (af::max<unsigned long>(maxs)>s.edge_select) break;
    }
    if (myWindow.close()) {
        cout<<endl<<"== Calculation stopped =="<<endl; 
        exit(1);
    }
    unsigned long mpos = posMax1D(maxs);
    myWindow(0,1).hist(af::constant(0.1f,maxs.dims(0)),0,maxs.dims(0));
    myWindow.show();
    af::deviceGC(); 
    return disks.slice(mpos);
}

float find_smooth(ifstream &fs, const DefaultVals &s, af::Window &myWindow, const af::array &disk0, float refine) {
    af::array maxs = af::constant(0.0,s.number_smooths);
    af::array disk, d_edge, disk0_smooth, edge0_smooth, ccor, smooths;
    while (!myWindow.close()) {
        disk = get_disks_RANDOM(fs,s,1);

        disk = af::medfilt(disk);
        d_edge = generate_edge(disk,1.0);
        d_edge = af::medfilt(d_edge);
        if (s.hann) d_edge = d_edge * generate_hann(d_edge);
        d_edge = normalizeCC(d_edge);
                
        float locmax;
        float globmax = 0;
        unsigned int pos;
        //this could be vectorised but not so important
        //for (int x=0;x<maxs.dims(0);x++) {
    
        smooths = af::range(s.number_smooths);
        if (refine < 1e-6) {
            smooths = 0.2f+smooths/2.0f;///5.0f;
        } else {
            smooths = 2 * refine / s.number_smooths * smooths;
        }
        disk0_smooth = smooth(disk0,1,smooths);
        edge0_smooth = generate_edge(disk0_smooth,1.0);
        edge0_smooth = af::medfilt(edge0_smooth);
        if (s.hann) edge0_smooth = edge0_smooth * generate_hann(edge0_smooth);
        edge0_smooth = normalizeCC(edge0_smooth);
        af::array d_edge_tiled = af::tile(d_edge,1,1,s.number_smooths);
        ccor = cross_correlate2(d_edge_tiled,edge0_smooth);
        int max_pos = posMax1D(af::flat(af::max(af::max(ccor,0),1)));
        maxs(seq(max_pos,max_pos)) =  maxs(seq(max_pos,max_pos)) + 1.0;
                
        myWindow(0,0).image(normaliseView(d_edge));
        myWindow(0,1).hist(100*maxs,0,maxs.dims(0),"smooth value");
        myWindow(0,2).image(normaliseView(smooth(disk0,0.2f+af::max<unsigned long>(maxs)/5.0f)),"tested smooth");

        myWindow.show();

        if (af::max<unsigned long>(maxs)>s.smooth_select) break;
    }
    if (myWindow.close()) {
        cout<<"== Calculation stopped =="<<endl; 
        exit(0);;
    }
    unsigned int pos_smooth = posMax1D(maxs);
    af::deviceGC();
    return af::sum<float>(smooths(seq(pos_smooth,pos_smooth)));
}

af::array finalise_disk(const DefaultVals &s, const af::array &disk0) {
    //upscale disk, smooth remove zagging edge smooth with found smooth0 and downscale for DPC processing
    af::array disk0_up = af::resize(disk0,5*disk0.dims(0),5*disk0.dims(1),AF_INTERP_BILINEAR);
    float sumdisk = af::sum<float>(disk0_up);
    disk0_up = smooth(disk0_up,10.0);
    disk0_up = maxnorm(disk0_up);
    float thr = 0.0f;
    while (true) {
        float current = af::sum<float>(1.0f*(disk0_up>thr));
        if (current < sumdisk) break;
        thr = thr + 0.01;
    }
    disk0_up = 1.0f*(disk0_up>thr);
    disk0_up = smooth(disk0_up,5*s.smooth_precise);

    af::array disk_final = af::resize(disk0_up, s.SIZE_DET * s.DOWNSCALE, s.SIZE_DET * s.DOWNSCALE,AF_INTERP_BILINEAR);
    disk_final = smooth(disk_final,s.smooth_precise);
    disk_final = maxnorm(disk_final);

    //if (s.ZEROPAD > 1) return zeropad(s,disk_final);
    //else 
    return disk_final.as(f32);
}

af::array removeHotPixels(const af::array &in, float thr) {
    //compute histogram
    af::array hist = af::histogram(in,256);
    // remove hotest pixels 0.01% -ish
    // CUDA code faills if the following or a hanning window is not used
    float sumhist = af::sum<float>(hist);
    int r = hist.dims(0)-1;
    af::array subhist;
    //get rid of 0.1 % of hot pixels
    while (2>1) {
        subhist = hist(seq(r--,hist.dims(0)-1));
        float sumsub = af::sum<float>(subhist);
        if (sumsub/sumhist > thr) break;
    }
    
    float min_in = af::min<float>(in);
    float max_in = af::max<float>(in);
    float threshold = (max_in - min_in)*r/hist.dims(0) + min_in;
    
    af::array out = in.copy();
    out = (in<threshold)*in + (in>=threshold)*threshold;
    af::deviceGC();
    return out;
}

void plot1d(af::Window &win, const af::array &ar1d, std::string lala) {
    //plot 1D line graph
    af::array x = af::range(ar1d.dims(0));
    win.plot(x.as(f32),ar1d.as(f32),lala.c_str()); //scatter with AF_MARKER_CIRCLE
}

af::array cross_correlate2(const af::array &a1, const af::array &a2) { 
    //cross correlation with upscale
    unsigned int scale = 2;
    af::array A1 = af::fft2(a1,a1.dims(0)*scale,a1.dims(1)*scale);
    af::array A2 = af::fft2(a2,a2.dims(0)*scale,a2.dims(1)*scale);
    af::array ccor_line = af::real(af::ifft2(A1*af::conjg(A2)));
	// INVERT QUADRANTS AFTER IFFT
	ccor_line = shift(ccor_line,ccor_line.dims(0)/2,ccor_line.dims(1)/2);
    ccor_line /= float(ccor_line.dims(0))*float(ccor_line.dims(1));
    unsigned int l = a1.dims(0)/2*(scale - 1);
    unsigned int h = a1.dims(0)/2*(scale + 1) -1;
    return ccor_line(seq(l,h),seq(l,h));
}


af::array get_disks_RANDOM(std::ifstream& fs, DefaultVals s,  unsigned long noFrames) {
    //calculate size of read out data for the *tmp0 pointer in bytes

    int FRAMELENGTH = s.HEADERSIZE+s.BYTESIZE_DET+s.FOOTERSIZE;
    int LENGTH = FRAMELENGTH * noFrames; 
    int DATALENGTH = s.BYTESIZE_DET * noFrames;

    //stack overflow answer to random number generation
    std::mt19937 rng;
    rng.seed(std::random_device()());
    // distribution to search for random set of frames not over the edge of acquisition
    unsigned long limit = s.SCAN_X;
    if (limit>s.SCAN_Y) limit = s.SCAN_Y;
    std::uniform_int_distribution<std::mt19937::result_type> dist6(0,limit - noFrames); 
    
    unsigned long rx = dist6(rng);
    unsigned long ry = dist6(rng);
    //somehow the random number overflows 
    while (rx>limit-noFrames) rx = dist6(rng);
    while (ry>limit-noFrames) ry = dist6(rng);
    
    unsigned long newpos = FRAMELENGTH * rx + FRAMELENGTH * ry * s.SCAN_X;

    fs.seekg(newpos, ios::beg);
    
    char *tmp0 = new char [LENGTH];
    fs.read((char*)tmp0, LENGTH);

    // 16-bit or 32-bit endian correction
    if (s.endian) endianSwap(s, LENGTH, tmp0);
  
    char *tmpNoHead = new char [DATALENGTH];
    //get rid of headers and footers by memcpy 
    for (int j = 0; j<noFrames;j++) {
        memcpy(tmpNoHead + j*s.BYTESIZE_DET, tmp0 + s.HEADERSIZE + j*FRAMELENGTH, s.BYTESIZE_DET);
    }

  	//read whole line into a new 3D array
  	s.STEP = noFrames;
   	af::array ar_line_small = pointerToArray(s, tmpNoHead);

    delete[] tmp0;
    delete[] tmpNoHead;     
     
    //ar_line_small = ar_line_small * af::tile(generate_hann(ar_line_small),1,1,ar_line_small.dims(2));

    //we can do this as averaging trick (smear in x direction)
    //ar_line_small = ar_line_small + af::shift(ar_line_small,0,0,1)/2.0 + af::shift(ar_line_small,0,0,-1)/2.0;
    if (s.DOWNSCALE < 0.99) ar_line_small = downscale(s,ar_line_small);
    //if (s.ZEROPAD > 1) ar_line_small = zeropad(s,ar_line_small);
    return ar_line_small;
}

std::string sConvert (float number)
{
    std::ostringstream buff;
    buff<<number;
    return buff.str();
}

void endianSwap(const DefaultVals &s, int LENGTH, char* tmp) {
    //char*c0 = (char*)tmp0;
    //for (int i = 0; i < LENGTH*2; i += 2)
    //swap(c0[i], c0[i + 1]);
    // using GCC builtin swap because it is about 30% faster 
    if (s.bytedepth == 2) {
        unsigned short*c0 = (unsigned short*)tmp;
        for (int i = 0; i < LENGTH/2; i ++)
            c0[i] = __builtin_bswap16(c0[i]);
    }
    if (s.bytedepth == 4) {
        float*c0 = (float*)tmp;
        for (int i = 0; i < LENGTH/4; i ++)
            c0[i] = __builtin_bswap32(c0[i]);
    }
}

af::array pointerToArray(DefaultVals s, char* tmpNoHead, bool mask) {
    af::array ar_line_small;
    if (mask) s.STEP = 1;
    if (s.bytedepth == 1) {
        if (s.sign) ar_line_small = af::array(s.SIZE_DET,s.SIZE_DET,s.STEP, (char*) tmpNoHead);
        else ar_line_small= af::array(s.SIZE_DET,s.SIZE_DET,s.STEP, (unsigned char*) tmpNoHead);        
    }    
    if (s.bytedepth == 2) {
        if (s.sign) ar_line_small = af::array(s.SIZE_DET,s.SIZE_DET,s.STEP, (short*) tmpNoHead);
        else ar_line_small =  af::array(s.SIZE_DET,s.SIZE_DET,s.STEP, (unsigned short*) tmpNoHead);         
    }
    if (s.bytedepth == 4) {
        ar_line_small= af::array(s.SIZE_DET,s.SIZE_DET,s.STEP, (float*) tmpNoHead);
    }
    return ar_line_small;
}

af::array downscale(const DefaultVals &s,const af::array &ar_line_small) {
    return af::resize(ar_line_small,int(ar_line_small.dims(0)*s.DOWNSCALE),int(ar_line_small.dims(1)*s.DOWNSCALE),AF_INTERP_BILINEAR);
}

af::array zeropad(const DefaultVals &s, const af::array &ar_line_small) {
    unsigned long newSize = ar_line_small.dims(0) * s.ZEROPAD; 
    unsigned long L = (newSize - ar_line_small.dims(0)) / 2;
    unsigned long H = L + ar_line_small.dims(0) - 1 ;
    af::array ar_line(newSize,newSize,ar_line_small.dims(2),f32);
    ar_line=0.0;
    ar_line(seq(L,H),seq(L,H),seq())=ar_line_small.as(f32);
    return ar_line;
}

array readDiskLine(ifstream& fs, const DefaultVals &s) {
    //calculate size of read out data for the *tmp0 pointer in bytes

    int FRAMELENGTH = s.HEADERSIZE+s.BYTESIZE_DET+s.FOOTERSIZE;
    int LENGTH = FRAMELENGTH * s.STEP; 
    int DATALENGTH = s.BYTESIZE_DET * s.STEP;

    char *tmp0 = new char [LENGTH];
    fs.read((char*)tmp0, LENGTH);

    
    //printf( "Here: %p\n" , tmp0 );

    // 16-bit or 32-bit endian correction
    if (s.endian) endianSwap(s, LENGTH, tmp0);
  
    char *tmpNoHead = new char [DATALENGTH];
    //get rid of headers and footers by memcpy 
    for (int j = 0; j<s.STEP;j++) {
        memcpy(tmpNoHead + j*s.BYTESIZE_DET, tmp0 + s.HEADERSIZE + j*FRAMELENGTH, s.BYTESIZE_DET);
    }

  	//read whole line into a new 3D array
   	af::array ar_line_small = pointerToArray(s, tmpNoHead);

    delete[] tmp0;
    delete[] tmpNoHead;     
     
    //ar_line_small = ar_line_small * af::tile(generate_hann(ar_line_small),1,1,ar_line_small.dims(2));

    //we can do this as averaging trick (smear in x direction)
    //ar_line_small = ar_line_small + af::shift(ar_line_small,0,0,1)/2.0 + af::shift(ar_line_small,0,0,-1)/2.0;
    if (s.DOWNSCALE < 0.99) ar_line_small = downscale(s,ar_line_small);
    //if (s.ZEROPAD > 1) ar_line_small = zeropad(s,ar_line_small);
    return ar_line_small;
}

af::array generate_edge(DefaultVals s, af::array mask) {
    //edge creation
    af::array A,B,dA,dB,gKerr;
    if (s.DOG) {
        A = af::gaussianKernel( mask.dims(0), mask.dims(1),s.dogLarge,s.dogLarge);
        B = af::gaussianKernel( mask.dims(0), mask.dims(1),s.dogSmall,s.dogSmall);
        return af::fftConvolve2(A,mask)-af::fftConvolve2(B,mask);
    } else if (s.sigma > 0.19) {
        gKerr = af::gaussianKernel( mask.dims(0), mask.dims(1),s.sigma,s.sigma);
        grad(A,B,gKerr);     

        dA = fftConvolve2(A,mask,AF_CONV_DEFAULT );
        dB = fftConvolve2(B,mask,AF_CONV_DEFAULT );
        return af::hypot(dA,dB);
    } else {
        // most datasets do really well with simple gradient (it probably requires smoothed mask done by AUTO algorithm)
        grad(A,B,mask);     
        
        return af::hypot(A,B); // dogo filter can work but its not better dog(mask,4.0f,1.0f);
    }
}

af::array generate_edge(const af::array &disk, float smooth) {
     af::array gKerr,dx,dy;  
     unsigned int scale = 2;
     gKerr = af::gaussianKernel(disk.dims(0), disk.dims(1),smooth,smooth);
     grad(dx,dy,gKerr);     
     af::array mask_dx = fftConvolve2(dx,disk,AF_CONV_DEFAULT);
     af::array mask_dy = fftConvolve2(dy,disk,AF_CONV_DEFAULT);
     af::array edge = af::hypot(mask_dx,mask_dy);
     
     af::deviceGC();

     return edge;
}


af::array zoom_image(const af::array &disk, float zoom) {
    //create transformation matrix according to af documentation
    float tx = float(disk.dims(0))*(1.0-zoom)*0.5;
    float ty = float(disk.dims(1))*(1.0-zoom)*0.5;
    float hA[] = {zoom, 0.0, tx, 0.0, zoom, ty};
    af::array trans(3,2,hA); 
    //zoom the array
    af::array zoomed = af::transform(disk,trans,disk.dims(0),disk.dims(1),AF_INTERP_BILINEAR);
    return zoomed;
}

void gotHere() {
    cout << "\nGot here!!" << endl;
    exit(0);
}


void check_s(DefaultVals &s) {
    //dimension of computed arrays, number of parallel frames depends on this - ideally larger card will give result anyways

    if (s.DOG) if (!s.checkDOG1 || !s.checkDOG2) {
        cout << "\n both sigmas (-t -u) for DOG filter need to be specified\n\n";
        exit(1);
    }

    if (s.DOG) {
        if (s.dogSmall>s.dogLarge) {
            double dogBuff=s.dogSmall;
            s.dogSmall=s.dogLarge;
            s.dogLarge=dogBuff;
        }
        if (s.dogSmall==s.dogLarge) {
            cout << "Do no use the same sigma for DOG filtering - you get zero!\n\n"; 
            exit(1);
        }
    }

    s.SCAN_XY        = s.SCAN_X * s.SCAN_Y;
    s.bytedepth      = s.bitdepth / 8;
    s.BYTESIZE_DET   = s.bytedepth * s.SIZE_DET * s.SIZE_DET;
    //read data by for loop one s.STEP at a time
    s.leftovers = s.SCAN_XY % s.STEP;
}

af::array get_mask(const DefaultVals &s) {     
    af::array tmp = af::loadImage(s.maskfile.c_str());
    if (s.DOWNSCALE < 0.99) tmp = downscale(s,tmp);
    return tmp.as(f32);
}

void setDefaultSettings(DefaultVals &s) {
    //default values for analysis data needs to be 16bit unsigned - otherwise convert by imageJ or ..
    s.AUTO       = false;    // automatic edge/mask generation
    s.DOG        = false;	 // true for difference of gaussians to be used
    s.sigma      = 0.0;	   	 // sigma for dx,dy Gradient analysis
    s.dogSmall   = 0.50;	 // small sigma DOG filter
    s.dogLarge   = 0.10;	 // large sigma DOG filter
    s.STEP       = 256;		 // number of images analysed at same time (257 is prime number)
    s.HEADERSIZE = 384;	     // 384 new files and 256 old ones
    s.FOOTERSIZE = 0;	     // empad files have 1024 byte footers
    
    s.endian     = true;     // correct endians - if used in windows or data were acquired on different endian machine
    s.show       = true;     // show the disks - no speed difference with direct opengl connection
    
    s.dirname = "outputDir"; //dir_name and position of files
    
    //DEFAULT IMPORTANT DIMENSIONS OF SCAN ANALYSIS (adjust according the scan and needed 0 padding)
    s.ZEROPAD    = 2;       //doubles the size of the image by zeropadding to protect for close to edge effects in cross-correlation processing and Nyqist errors
    s.DOWNSCALE  = 1.0;     // for faster but less precise analysis
    s.SIZE_DET   = 256;     //256 for Medipix3 single
    s.SCAN_X     = 257;     //257
    s.SCAN_Y     = 256;     //256

    //select if hanning (not needed by default -> disk doesnt touch the edge!)
    //the implementation of 2D hanning is probably wrong
    s.hann = false;
    
    s.sign       = false;    //switch for signed/unsigned data from detector
    s.checkDOG1  = false;
    s.checkDOG2  = false;
    s.bitdepth   = 16;
    
    //AUTO SETTINGS
    s.number_edges = 64;
    s.edge_select = 32;
    s.number_smooths = 32;
    s.smooth_select = 32;
    s.summing = 32;
}


af::array partialConvolve(const DefaultVals &s, const af::array &edge_line) //done to speed up convolution of the edges - saving one fft which is done before sequential calculation
{   
    af::array Edge_line = fft2(edge_line,edge_line.dims(0)*s.ZEROPAD,edge_line.dims(1)*s.ZEROPAD);
    Edge_line = Edge_line*tile(s.Edge,1,1,Edge_line.dims(2));
	af::array ccor_line = real(ifft2(Edge_line));
	// INVERT QUADRANTS AFTER IFFT
	ccor_line = shift(ccor_line,ccor_line.dims(0)/2,ccor_line.dims(1)/2);
    unsigned int l = edge_line.dims(0)/2*(s.ZEROPAD - 1);
    unsigned int h = edge_line.dims(1)/2*(s.ZEROPAD + 1) -1;
    return ccor_line(seq(l,h),seq(l,h))/float(ccor_line.dims(0))/float(ccor_line.dims(1))*s.ZEROPAD*s.ZEROPAD;
}


af::array cross_correlate(const af::array &a1, const af::array &a2) 
{   
    af::array A1 = af::fft2(a1);
    af::array A2 = af::fft2(a2);
    af::array ccor_line = af::real(af::ifft2(A1*af::conjg(A2)));
	// INVERT QUADRANTS AFTER IFFT
	ccor_line = shift(ccor_line,ccor_line.dims(0)/2,ccor_line.dims(1)/2);

    return ccor_line/float(ccor_line.dims(0))/float(ccor_line.dims(1));
}


Results createResults(const DefaultVals &s) {
    //reinitiate results - if compute failed due to GPU memory overflow
    Results calcs;
    //generate arrays in 3D setting - easier to handle, will be reshaped later
    calcs.maxs = af::constant(0.0,s.SCAN_XY);
    calcs.sums = af::constant(0.0,s.SCAN_XY);
    calcs.poss = af::constant(0.0,s.SCAN_XY,2);
    calcs.subs = af::constant(0.0,s.SCAN_XY,4);
    calcs.dpcx = af::constant(0.0,s.SCAN_XY);
    calcs.dpcy = af::constant(0.0,s.SCAN_XY);    
    return calcs;
}
    
void generate_outputs(const DefaultVals &s,Results &calcs) {
    //sum and maximum of correlation image
    calcs.maxs = af::moddims(calcs.maxs, s.SCAN_X,s.SCAN_Y);
    calcs.sums = af::moddims(calcs.sums, s.SCAN_X,s.SCAN_Y);
    //DPC deflection images

    //subpixel calculation arrays
	array s01 = flat(calcs.subs(seq(),seq(0,0)));
	array s10 = flat(calcs.subs(seq(),seq(1,1)));
	array s11 = flat(calcs.maxs);
	array s12 = flat(calcs.subs(seq(),seq(2,2)));
	array s21 = flat(calcs.subs(seq(),seq(3,3)));

    //subpixel quadratic calculation
    calcs.dpcx = (s10-s12)/2.0/(s10+s12-2.0*s11);
    calcs.dpcy = (s01-s21)/2.0/(s01+s21-2.0*s11);	
    calcs.dpcx = af::moddims(calcs.dpcx, s.SCAN_X,s.SCAN_Y);
    calcs.dpcy = af::moddims(calcs.dpcy, s.SCAN_X,s.SCAN_Y);
	//adding integer pixel measure of deflection
	calcs.dpcx = calcs.dpcx + af::moddims(calcs.poss(seq(),seq(0,0)), s.SCAN_X,s.SCAN_Y);
	calcs.dpcy = calcs.dpcy + af::moddims(calcs.poss(seq(),seq(1,1)), s.SCAN_X,s.SCAN_Y);
}

void data_check(ifstream& fs, const DefaultVals &s) {
	//check for datafile existence
	if (!fs) {
	    cout << "\nSelected raw/bin/mib/* file does not exist!\n\n"; exit(0);
	}
	
    fs.seekg(0, ios::end);

	ulong filesize = fs.tellg();
    
    cout << "\nbytedepth = "<<s.bytedepth<<endl;
    cout << "header = " << s.HEADERSIZE << endl;
    cout << "footer = " << s.FOOTERSIZE << endl;
    cout << "detector = " << s.BYTESIZE_DET << endl;
    
    ulong dataframes = ulong(filesize/(s.HEADERSIZE+s.FOOTERSIZE+s.BYTESIZE_DET));
	cout << "\nSize of the file is: " << filesize << " B (" << filesize/1024/1024 << " MB / " << filesize/1024/1024/1024 <<"GB)";
	cout << "\nNumber of frames in binary file: " << dataframes;
	cout << "\nNumber of frames in xy parameters: " << s.SCAN_XY;
	  
    if (dataframes < s.SCAN_XY) {
	    cout << "\nFile is too small for the number of scan -x and -y parameters!\n\n";
	    exit(0);
	}
	
    //go back to begining
    fs.seekg(0, ios::beg);
}

void visualise_calc(DefaultVals &s, Results &results, af::Window& myWindow, const af::array &disk, const af::array &edge, const af::array &ccor, long indL, long indH) {
    af::array s01,s10,s11,s12,s21,a,e,f,g,h;
    //the live view needs to be normalised only up to current analysis point in 3D array
    //after that its only zeros
    s01 = flat(results.subs(seq(),seq(0,0)));
    s10 = flat(results.subs(seq(),seq(1,1)));
    s11 = flat(results.maxs);
    s12 = flat(results.subs(seq(),seq(2,2)));
    s21 = flat(results.subs(seq(),seq(3,3)));

    //subpixel quadratic calculation
    results.dpcx = (s10-s12)/2.0f/(s10+s12 - 2.0f * s11);
    results.dpcy = (s01-s21)/2.0f/(s01+s21 - 2.0f * s11);	
    //adding integer pixel measure of deflection
    results.dpcx = results.dpcx + results.poss(seq(),seq(0,0));
    results.dpcy = results.dpcy + results.poss(seq(),seq(1,1));
    
    a = results.sums.as(f32);
    a = normaliseViewRange1D(a,indH);
    a = af::moddims(a,s.SCAN_X,s.SCAN_Y);
    a = rotate(a,-M_PI/2.0f,false);
    
    e = results.dpcx.as(f32);
    e = normaliseViewRange1D(e,indH);
    e = af::moddims(e,s.SCAN_X,s.SCAN_Y);
    e = rotate(e,-M_PI/2.0f),false;

    f = results.dpcy.as(f32);
    f = normaliseViewRange1D(f,indH);
    f = af::moddims(f,s.SCAN_X,s.SCAN_Y);
    f = rotate(f,-M_PI/2.0f),false;
    
    g = results.maxs.as(f32);
    g = normaliseViewRange1D(g,indH);
    g = af::moddims(g,s.SCAN_X,s.SCAN_Y);
    g = rotate(g,-M_PI/2.0f),false;
    
    h = sqrt((e-0.5)*(e-0.5) + (f-0.5)*(f-0.5));
    h = normaliseViewRange1D(h,indH);
    h = af::moddims(h,s.SCAN_X,s.SCAN_Y);
    
    
    visualisation_window(a,
                normaliseView(disk),
                normaliseView(edge),
                normaliseView(ccor),
                e,f,g,h,
                myWindow);
}


af::array autoMask(ifstream& fs, DefaultVals &s) {
    timer::start();
    //arrays to store edges and disk
    af::array edges,disks;
    
    //create window
    af::Window WinAuto(3*256,32 + 256,"disk AUTO mask");
    WinAuto.grid(1,3);
    WinAuto.setSize(3*256,32 + 256);
    
    //function to generate disks and edges from random points of dataset
    generate_disks_edges(fs,s,WinAuto,disks,edges);
        
    //vote and select the best disk
    af::array disk0 = select_disk(fs,s,WinAuto,disks,edges);
    
    //look for a rough smooth value to match the data
    float smooth_rough = find_smooth(fs,s,WinAuto,disk0);
    //look for precise smooth value
    s.smooth_precise = find_smooth(fs,s,WinAuto,disk0, smooth_rough);
    
    cout << "\n\n\nStatistically evaluated smooth: sigma = " << s.smooth_precise<<endl;
    //generate the final mask - done in an upscaled version for accuracy
    af::array disk_final = finalise_disk(s,disk0);
    
    //tidy up
    WinAuto.setVisibility(false);
    af::sync();
    af::deviceGC();
    
    cout << "Mask autofind algorithm elapsed seconds: " << roundFloat(timer::stop()) << " s";
    
    //return mask variable
    return disk_final;
}

af::array mask_Edge(DefaultVals &s) {
    //generate edge
    s.edge = generate_edge(s,s.mask);
    //normalise for cross-correlation - this is actually not at all needed but doesnt cost much if done only once
    normalizeCCinPlace(s.edge);
    //prepare for edge phase cross-correlation - this will save some computation time 1 out of 5 ffts
	return conjg(fft2(s.edge,s.edge.dims(0)*s.ZEROPAD,s.edge.dims(1)*s.ZEROPAD));
}

int compute(ifstream& fs, DefaultVals &s, Results &results, af::Window& myWindow, bool end) {
    af::array ar_line, edge_line, diskmask, ccor_line, subs;
    float deltaT = 1.0;
    while (s.frame_i < s.SCAN_XY - s.leftovers || end) {  
        //this is sort of a funky way but end needs to be set true only once before the
        //calculation of leftovers and immediatelly switched back so while case stops
        end = false;

        //read a set of disks from the current position in fs stream
        ar_line = readDiskLine(fs,s); 

        //set low and high index for data storage navigation in 3D calc results
        long indL = s.frame_i;
        long indH = indL + s.STEP - 1;
       
        //sum data with a vector algorithm
        results.sums(seq(indL,indH)) = af::flat(af::sum(af::moddims(ar_line,ar_line.dims(0)*ar_line.dims(1),ar_line.dims(2)),0));

        //median filtering to remove any shot noise
        ar_line = af::medfilt(ar_line);
        
        //generate edges by DOG or grad Gaussian methods
        edge_line = generate_edge(s, ar_line.as(f32));

        //remove hot pixels at the edge of gradient
        diskmask = smoothed_disk(edge_line.dims(0),int(edge_line.dims(0)*0.45),5);
        edge_line = af::tile(diskmask,1,1,edge_line.dims(2)) * edge_line;
    
        //create cross-correlation pattern by phase convolution - s.Edge was pre-FFT-ed to save one FFT operation per cycle
        //edge_line = af::medfilt(edge_line);
        normalizeCCinPlace(edge_line);
        ccor_line = partialConvolve(s,edge_line);           

        af::array tmpmaxnew, tmpposnew;
        //find positions of max with a vector algorithm
        af::max(
            tmpmaxnew,tmpposnew,af::moddims(ccor_line,ccor_line.dims(0)*ccor_line.dims(1),ccor_line.dims(2)),0
        );
        results.maxs(seq(indL,indH)) = flat(tmpmaxnew);
        tmpposnew = flat(tmpposnew);
        
        //generate xy position of maxima
        results.poss(seq(indL,indH),seq(0,0)) = af::mod(tmpposnew,ccor_line.dims(0));
        results.poss(seq(indL,indH),seq(1,1)) = tmpposnew / ccor_line.dims(0);
        
        //read out surrounding pixel for subpixel registration
        //this is done by a lookup code - we know positions of max so we can look up one to the left,top,bottom and right
        //this way the code finds position of each maxima globaly
        tmpposnew = tmpposnew + seq(0,ccor_line.dims(0)*ccor_line.dims(1)*s.STEP-1,ccor_line.dims(0)*ccor_line.dims(1));
        //this array contains indices of surrounding pixels
        tmpposnew = af::join(1,tmpposnew-ccor_line.dims(0),tmpposnew-1,tmpposnew+1,tmpposnew+ccor_line.dims(0));
        //lookup the values of surrounding pixels
        subs =  af::lookup(af::flat(ccor_line),flat(reorder(tmpposnew,1,0)));
        //reorder them in a right way - this was trial and error -> things get complicated
        results.subs(seq(indL,indH),seq()) = af::reorder(af::moddims(subs,4,s.STEP),1,0);    

        //show the calculation 
        if (s.show) {
            visualise_calc(s, results, myWindow, ar_line.slice(0), edge_line.slice(0), ccor_line.slice(0), indL, indH);
        }
        
        s.frame_i+=s.STEP;
        
        // do some estimation for finishing time (hopefully better than windows)
        float perc = 100.0*(float(s.frame_i)/float(s.SCAN_XY));
        float dt = timer::stop(s.start);
        float framerate = s.STEP /  (dt - deltaT);
        deltaT = dt;
        dt = dt / perc * 100.0 - dt;
        cout << "\r"<< roundFloat(perc) << "%      \t finishing in ~ "  << roundFloat(dt) << " s     \t" << roundFloat(framerate) << " fps      " << s.frame_i << " done out of " << s.SCAN_XY << "        ";
        
        //tidy up
        af::deviceGC();
    if (s.show and myWindow.close()) {
        cout << "\r=== CALCULATION TERMINATED === " << endl; 
        exit(1);
    }
    //allow for termination by closing the window - if opened
    }
}


void visualisation_window(const af::array &a, const af::array &b, const af::array &c, const af::array &d, const af::array &e, const af::array &f, const af::array &g, const af::array &h, af::Window& myWindow) {
    myWindow.grid(2,4);
    myWindow(0,0).image(a, "sum");
    myWindow(0,1).image(b, "raw disk");
    myWindow(0,2).image(c, "edge");
    myWindow(0,3).image(d, "cross-corr");    
    myWindow(1,0).image(e, "dpcX");
    myWindow(1,1).image(f, "dpcY");
    myWindow(1,2).image(g, "corr-val");
    myWindow(1,3).image(h, "magnitude");

    myWindow.show();
}


void createLogFile(const DefaultVals &s) {
    //GENERATE LOGFILE RECORDING USED PARAMETERS
	string logfile="./log.txt";
    ofstream outputFile(logfile.c_str()); //this may not require c string conversion

    outputFile << "_____________________________________" << endl;
    outputFile << "__________ANALYSIS LOG FILE__________" << endl;
    outputFile << "_____________________________________" << endl;

    outputFile  <<"binary file: = -b "<< s.binfile <<endl;
    outputFile  <<"output directory: -o "<<s.dirname <<endl;
    if (s.AUTO) {
        outputFile << "automatic mask generation used -A"<<endl;
        outputFile << "\tmask saved as mask.tiff"<<endl;
        outputFile << "\tnumber of frames averaged for analysis = "<<s.summing<<endl;
        outputFile << "\tnumber of edges generated for statistical comparison = "<<s.number_edges<<endl;
        outputFile << "\tnumber of votes to select an edge  = "<<s.edge_select<<endl;
        outputFile << "\tnumber of different smooth settings evenly spread to explore = "<<s.number_smooths<<endl;
        outputFile << "\tnumber of votes to select a smooth value = "<<s.smooth_select<<endl;
        outputFile << "\tthresholded mask smothed with sigma = "<<s.smooth_precise<<endl;
    } else {
        outputFile  <<"mask file: "<< s.maskfile <<endl;
    }

    if (s.DOG) {
        outputFile  <<"difference of Gaussian filter used"<<endl;
        outputFile  <<"sigma1: "<< s.dogSmall <<endl;
        outputFile  <<"sigma2: "<< s.dogLarge <<endl;
    }
    else {
            outputFile  <<"Gaussian edge convolution filter used"<<endl;
            outputFile  <<"sigma: "<< s.sigma <<endl;
    }

    outputFile<<"output dir: "<< s.dirname <<endl; //this is fairly obvious
    outputFile<<"header size: "<< s.HEADERSIZE <<" kb" <<endl;
    outputFile<<"footer size: "<< s.FOOTERSIZE <<" kb" <<endl;
    outputFile<<"number of files computed in parallel: "<< s.STEP <<endl;
    outputFile<<"scan dimensions: X = "<< s.SCAN_X<<"  Y = "<< s.SCAN_Y <<endl;
    outputFile<<"size of detector (square): "<< s.SIZE_DET <<endl;
    outputFile<<"zero padding: "<< s.ZEROPAD <<endl;
    outputFile<<"downscale: "<< s.DOWNSCALE <<endl;
    outputFile<<"calculation time (seconds): "<< timer::stop() <<endl;
    outputFile.close();
}

int getSettings(DefaultVals &s, int argc, char** argv) {
    //this looks bad but it is simple
    while ((s.ind = getopt (argc, argv, "S:N:M:c:O:T:m:e:b:p:o:x:y:z:a:d:s:t:u:hwlgfAF")) != -1)
        switch (s.ind) {
            case 'h':
                showhelpinfo();
                exit(1);
            case 'm':
                s.maskfile = optarg;
                break;
            case 'e':
                s.HEADERSIZE = atoi(optarg);
                break;
            case 'p':
                s.STEP = atoi(optarg);
                break;
            case 'd':
                s.SIZE_DET = atoi(optarg);
                break;
            case 'x':
                s.SCAN_X = atoi(optarg);
                break;
            case 'y':
                s.SCAN_Y = atoi(optarg);
                break;
            case 'a':
                s.bitdepth = atoi(optarg);
                break;
            case 'z':
                s.ZEROPAD = atoi(optarg);
                break;
            case 'b':
                s.binfile = optarg;
                break;
            case 'o':
                s.dirname = optarg;
                break;        
            case 'w':
                s.endian = false;
                break;
            case 'A':
                s.AUTO = true;
                break;
            case 'g':
                s.DOG=true;
                break;
            case 'f':
                s.FOOTERSIZE = atoi(optarg);
                break;
            case 'l':
                s.sign=true; 
                break;            
            case 's':
                s.sigma =  atof(optarg);
                break;            
            case 'c':
                s.DOWNSCALE =  1.0f / atof(optarg);
                break;
            case 't':
                s.tmp11 = optarg;
                s.dogLarge =  atof(s.tmp11.c_str());
                s.checkDOG1 = true;
                break;
            case 'u':
                s.tmp22 = optarg;
                s.dogSmall =  atof(s.tmp22.c_str());
                s.checkDOG2 = true;
                break;
            case 'S':
                s.summing = atoi(optarg);
                break;
            case 'N':
                s.number_edges = atoi(optarg);
                break;
            case 'M':
                s.edge_select = atoi(optarg);
                break;
            case 'O':
                s.number_smooths = atoi(optarg);
                break;
            case 'T':
                s.smooth_select = atoi(optarg);
                break;
            case 'F':
                s.number_edges = 8;
                s.edge_select = 16;
                s.number_smooths = 16;
                s.smooth_select = 4;
                s.summing = 8;
                break;
            default:
                showhelpinfo();
                return 1;
        }
}

void saveimages(const DefaultVals &s, Results calcs) {  
    //BUILD RESULTING IMAGES
    generate_outputs(s,calcs);
    
    //save generated edge_array for a reference
    saveImageNative("./edge.tif",s.edge);
    saveImageNative("./mask.tif",s.mask);
    
    //rotate the arrays - to correct detector to c++ scan rotation
    // shifts are changed due to the new pixel size if downscaled
    calcs.dpcx=rotate(calcs.dpcx,-M_PI/2.0f,false)/s.DOWNSCALE; 
    calcs.dpcy=rotate(calcs.dpcy,-M_PI/2.0f,false)/s.DOWNSCALE;
    calcs.maxs=rotate(calcs.maxs,-M_PI/2.0f,false);
    calcs.sums=rotate(calcs.sums,-M_PI/2.0f,false);

    //write the arrays

    //32 bit results
    saveImageNative("./ccor_max.tif",calcs.maxs);
	saveImageNative("./x.tif",calcs.dpcy);
	saveImageNative("./y.tif",calcs.dpcx);
	saveImageNative("./sum.tif",calcs.sums);

    calcs.dpcx-=mean<float>(calcs.dpcx);
    calcs.dpcy-=mean<float>(calcs.dpcy);
    
    //16bit results
	calcs.maxs=  normalizeImage(calcs.maxs);
	calcs.dpcx = normalizeImage(calcs.dpcx);
	calcs.dpcy = normalizeImage(calcs.dpcy);
	calcs.sums=  normalizeImage(calcs.sums);

	saveImageNative("./16bit_ccor_max.tif",calcs.maxs.as(u16));
	saveImageNative("./16bit_x.tif",calcs.dpcx.as(u16));
	saveImageNative("./16bit_y.tif",calcs.dpcy.as(u16));
	saveImageNative("./16bit_sum.tif",calcs.sums.as(u16));

    createLogFile(s);
      
    mkdir("pngs", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    chdir("pngs");
  
    // maximise contrast by histogram calculation and putting % pixel out of resulting png and shown images
    float percentage = 1.0;
    calcs.maxs = norm_hist(calcs.maxs,percentage);
    calcs.dpcx = norm_hist(calcs.dpcx,percentage);
    calcs.dpcy = norm_hist(calcs.dpcy,percentage);
    calcs.sums = norm_hist(calcs.sums,percentage);
    
    normalizeImageInPlace(calcs.maxs);
	normalizeImageInPlace(calcs.dpcx);
	normalizeImageInPlace(calcs.dpcy);
	normalizeImageInPlace(calcs.sums);
  
    saveImageNative("./16bit_ccor_max.png",calcs.maxs.as(u16));
	saveImageNative("./16bit_x.png",calcs.dpcx.as(u16));
	saveImageNative("./16bit_y.png",calcs.dpcy.as(u16));
	saveImageNative("./16bit_sum.png",calcs.sums.as(u16));
    
	//CREATING PREVIEW AND SHOWING IT STRAIGHT AFTER FINISHED
	array fehViewT = join(1,calcs.maxs.as(u16),calcs.sums.as(u16));
	array fehViewB = join(1,calcs.dpcx.as(u16),calcs.dpcy.as(u16));
	array fehView  = join(0,fehViewT,fehViewB);
    
	saveImageNative("./fehPreView.tif",fehView);
	system("feh ./fehPreView.tif &");
    
    chdir("..");
    cout << "Data saved to DIR = " << s.dirname << endl << endl;
    calcFinitoString(); 
}


void showCreator() {
   /* 
    printf(R"EOF(
██████╗ ██╗██╗  ██╗██████╗ ██████╗  ██████╗    ███████╗████████╗███████╗███╗   ███╗
██╔══██╗██║╚██╗██╔╝██╔══██╗██╔══██╗██╔════╝    ██╔════╝╚══██╔══╝██╔════╝████╗ ████║
██████╔╝██║ ╚███╔╝ ██║  ██║██████╔╝██║         ███████╗   ██║   █████╗  ██╔████╔██║
██╔═══╝ ██║ ██╔██╗ ██║  ██║██╔═══╝ ██║         ╚════██║   ██║   ██╔══╝  ██║╚██╔╝██║
██║     ██║██╔╝ ██╗██████╔╝██║     ╚██████╗    ███████║   ██║   ███████╗██║ ╚═╝ ██║
╚═╝     ╚═╝╚═╝  ╚═╝╚═════╝ ╚═╝      ╚═════╝    ╚══════╝   ╚═╝   ╚══════╝╚═╝     ╚═╝
)EOF");*/
	cout << "\nPixelated DPC\n";
	cout << "GPU accelerated software to analyse STEM disk shifts\n(C) Matus Krajnak 2018\n\n";
}

void calcString() {
    cout << "\n\n=========================\n";
    cout <<     "== CALCULATION STARTED ==\n";
    cout <<     "=========================\n\n";
}

void calcFinitoString() {
    cout << "==========================\n";
    cout << "== CALCULATION FINISHED ==\n";
    cout << "==========================\n\n";
}
