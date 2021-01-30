// imagen.cpp
// g++ imagen.cpp -o imagen -fopenmp && ./imagen 
// Reading a PPM file a convert to PGM file
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <string.h>

using namespace std;

typedef struct{
  int Red;
  int Green;
  int Blue;
} Pixel;

class Image{
  public:
  Pixel *RGB;
  int height, width;
  string name, type, comment, max_pixelvalue;
  // Constructors
  // 1: Read from PPM
  Image(char*);
  // Destructor
  ~Image();
  // Methods
  int *NewMatrix(int, int);
  void Zeros(int **, int);
  Pixel *NewPixel(int, int);
  void ZerosRGB(Pixel **, int);
  void Read_PPMImage(char*); 
  void PPMtoPGM();
  void GrayScale(int **);
  void ExportPGM(int *);
};

int main () {
/* clock_t clock(void) returns the number of clock ticks 
       elapsed since the program was launched.To get the number  
       of seconds used by the CPU, you will need to divide by  
       CLOCKS_PER_SEC.where CLOCKS_PER_SEC is 1000000 on typical 
       32 bit system.  */
    char handle[] = "luna239092020.ppm"; 
    clock_t start, end;
    /* Recording the starting clock tick.*/
    start = clock(); 
    // Task
    Image im1(handle);
    im1.PPMtoPGM();  
    // Recording the end clock tick. 
    end = clock(); 
    // Calculating total time taken by the program. 
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    cout << "Time : " << fixed  << time_taken << setprecision(5); 
    cout << " sec " << endl; 
  return 0;
}

// Constructor
Image::Image(char *handle){
  Read_PPMImage(handle);
}
// Destructor
Image::~Image(){
  // Free memory allocation
  delete []RGB;
}

// New instance for an array
int *Image::NewMatrix(int height, int width){
  int *M = new int [height*width];
  Zeros(&M, height*width);
  return M;
}
// Fill all int array with zeros
void Image::Zeros(int **M, int size){
  for(int i=0; i<size; i++){
    (*M)[i] = 0;
  }
}
// New instance for an pixel array
Pixel *Image::NewPixel(int height, int width){
  Pixel *M = new Pixel [height*width];
  ZerosRGB(&M, height*width);
  return M;
}
// Fill al rgb values pixel array with zeros
void Image::ZerosRGB(Pixel **M, int size){
  for(int i=0; i<size; i++){
    (*M)[i].Red = 0;
    (*M)[i].Green = 0;
    (*M)[i].Blue = 0;
  }
}
// Read a PPM file
void Image::Read_PPMImage(char *handle){
  string line;
  char *s;
  // control pixels
  int index = 0;
  // index matrix
  int i=0;
  int option = 1;
  int value;
  // Save the name file
  name += handle;

  ifstream myfile(handle);
  if (myfile.is_open()){
    // Type: P3
    getline(myfile, type);
    // Comment
    getline(myfile, comment);
    // Width and Height
    getline(myfile,line);
    sscanf (line.c_str(),"%d %d", &width, &height);
    // Max pixel value
    getline(myfile, max_pixelvalue);
    // Set dinamic memory
    RGB = NewPixel(height, width);

    while(getline(myfile,line)){
      // count three lines and change matrix index
      if(index%3==0){
        i++;
      }
      sscanf(line.c_str(),"%d", &value);
      // Set pixel in respective matrix value
      switch(option){
        case 1: // Red pixel
          RGB[i].Red = value;
          option = 2;
          break;
        case 2:  // Green pixel
          RGB[i].Green = value;
          option = 3;
          break;
        case 3: // Blue pixel
          RGB[i].Blue = value;
          option = 1;
          break;
        default: // delfault
          value = 0;
      }
      index++;
    }
  }   else cout << "Unable to open file";
}
// Convert PPM to PGM
void Image::PPMtoPGM(){
  // Declare new int array
  int *Gray = NewMatrix(height, width);
  // Calculate gray scale
  GrayScale(&Gray);
  // Write PGM file
  ExportPGM(Gray);
  // Free memory allocation
  delete []Gray;
}
// Calculate parallel gray scale transform
void Image::GrayScale(int **Gray){
  #pragma omp parallel for
  for(int i=0; i<height*width; i++){
    (*Gray)[i] = 0.299*RGB[i].Red + 0.587*RGB[i].Green + 0.114*RGB[i].Blue;
  }
}
// Write a PGM file
void Image::ExportPGM(int *Gray){
  // Set name of pgm file
  string pgmfile = name.substr(0, name.find(".")) + ".pgm";
  // write pgm file
  ofstream myfile (pgmfile.c_str());
  if (myfile.is_open()){
    // Header
    myfile << "P2\n";
    // Comment
    myfile << comment << "\n";
    // Image dimension 
    myfile << width << " " << height << "\n";
    // Max Pixel Value
    myfile << max_pixelvalue << "\n";
    // Set all values
    for(int i=0; i<height*width; i++){
      myfile << Gray[i] << "\n";
    }
    myfile.close();
  }
  else cout << "Unable to open file";
}