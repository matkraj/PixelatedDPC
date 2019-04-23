# Pixelated DPC 

GPU accelerated analysis of disk image position in pixelated STEM data with subpixel resolution. 
This code was developped for magnetic imaging in scanning transmission electron microscopy with fast pixelated detectors. It is based on a gradient cross-correlation with a mask pattern. The mask can be created automatically from the dataset or input  manually by a user. It was tested with datasets larger than 100GB. Any modern gaming computer should be suitable. The research is described in [paper (free access)](https://doi.org/10.1016/j.ultramic.2016.03.006) and  [thesis](http://theses.gla.ac.uk/7906/). If you use Pixelated DPC in your research, we kindly ask that you cite the following paper: 
> Krajnak, Matus, et al. "Pixelated detectors and improved efficiency for magnetic imaging in STEM differential phase contrast." Ultramicroscopy 165 (2016): 42-50.

(tba: step by step guide)

##### Requirements: 

1. Linux operating system with modern GPU and CUDA/OpenCL drivers, build system and cmake
2. Arrayfire library https://arrayfire.com/ (installation guide http://arrayfire.org/docs/installing.htm)

##### Pixelated DPC Installation:

```bash
$ git clone https://github.com/matkraj/pixelatedDPC.git
$ cd pixelatedDPC/src/
$ mkdir build && cd build
$ cmake ..
$ make
```
##### Usage:
Help file describing all analysis settings:
```bash
$ ./pixDPC -h
```
Simple example
```bash
$ ./pixDPC -b binaryFILE -o outputDIR -d 256 -x 512 -y 512 -e 384 -f 0 -a 16 -p 64 -A
```
> This will analyse a dataset from detector of size 256x256 which has a header information of 384 bytes in front of each frame. The scan is 512x512, bitdepth of the detector is 16 bit and the mask will be found by automatic statistical analysis of the data. The data will be analysed with 64 frames in parallel.

##### Speed tests (GTX 1080 Ti analysis framerates for 256x256 16bit detector):
- upscale 2x ~ 2300fps
- no change ~ 4700fps
- 2x binning ~ 7300fps

##### Notes:

- Code supports 8/16/32-bit binary streams. 
- Tested with datasets from Medipix3 (MerlinEM) and EMPAD detectors
- Switch for signed and unsigned images and a change of endian architecture 
