//============================================================================
// Name        : Dip1.cpp
// Author      : Ronny Haensch
// Version     : 2.0
// Copyright   : -
// Description : 
//============================================================================

#include "Dip1.h"

//funktion that sets a color for a range of pixels
/*
image image that schould be changed Mat&
point where the range of pixels begin Point
size the width and height of  the pixel-range int
color the color which should be set for the pixel-range Vec3b
*/
void setAvgColor(Mat& image, Point point, int size, Vec3b color) {
	int anzColumns = image.cols;
	int anzRows = image.rows;
	int endOfX = point.x + size;
	int endOfY = point.y + size;
	//walk through every column
	for (int i = point.x; i < endOfX; i++) {
		//proof for image bound in x
		if (i < anzColumns) {
			//walk through every row per column
			for (int j = point.y; j < endOfY; j++) {
				//proof for image bound in y
				if (j < anzRows) {
					//set everage colour
					image.at<Vec3b>(Point(i,j)) = color;
				}
			}
		}
	}
	return;
}

//funktion that pixelate an area 
/*
sp startPixel Point
size side length of one output pixel
return color of output pixel
*/
void pixelate(Mat& inputImage, Mat& outputImage, int sizeOfPixel) {
	int anzRows = inputImage.rows; //#rows in matrix
	int anzCols = inputImage.cols; //#columns in matrix
	Point point; //begin of a Pixel-range in outputImage
	//int sizeOfPixel = 2; //size of a pixel-range in outputImage
	Vec3b color;

	//walk through every column
	for (int i = 0; i < anzCols; i++) {
		//walk through every column
		if (i % sizeOfPixel == 0) {
			//walk through every row
			for (int j = 0; j < anzRows; j++) {
				if (j % sizeOfPixel == 0) {
					point.x = i;
					point.y = j;
					//color = getAvgColor(inputImage, point, sizeOfPixel);
					//color set by first pixel of group
					color = inputImage.at<Vec3b>(point);
					//set color for a group of pixels
					setAvgColor(outputImage, point, sizeOfPixel, color);
				}
			}
		}
	}
}

//funktion that change the color for every pixel of an image-Matrix
//the change depends on the spezific previous color of every pixel
/*
image picture that should be changed
newColor color that should be added for every pixel
*/
void stainPicture(Mat& image, double newColor[3]) {
	int anzRows = image.rows; //#rows in matrix
	int anzCols = image.cols; //#columns in matrix
	double oldColor[] = {0.0, 0.0, 0.0};
	int colorA;
	int colorB;

	cout << "BGR:" << newColor[0] << "," << newColor[1] << "," << newColor[2] << endl;
	//walk through every column
	for (int i = 0; i < anzCols; i++) {
		//walk through every row
		for (int j = 0; j < anzRows; j++) {
			//obtain the color blue/green/red from the current pixel
			oldColor[0] = image.at<Vec3b>(Point(i, j)).val[0];	//blue
			oldColor[1] = image.at<Vec3b>(Point(i, j)).val[1];	//green
			oldColor[2] = image.at<Vec3b>(Point(i, j)).val[2];	//red
			
			//set the new color {yellow, magenta, cyan} to the current pixel
			for (int k = 0; k < 3; k++){
				//proof for bounds of color-space
				//if beneath bounds -> 0 min
				if ((oldColor[k] - newColor[k]) < 0){
					image.at<Vec3b>(Point(i, j)).val[k] = 0;
				}
				//if above bounds -> 255 max
				else if ((oldColor[k] - newColor[k]) > 255) {
					image.at<Vec3b>(Point(i, j)).val[k] = 255;
				}
				//if not out of bounds -> calculate
				else {
					image.at<Vec3b>(Point(i, j)).val[k] = oldColor[k] - newColor[k];
				}
			}
		}
	}
	//image.at<Vec3b>(Point(i, j)).val[0] = oldColor[0] - newColor[0];
	//image.at<Vec3b>(Point(i, j)).val[1] = oldColor[1] - newColor[1];
	//image.at<Vec3b>(Point(i, j)).val[2] = oldColor[2] - newColor[2];
	return;
}

// function that performs some kind of (simple) image processing
/*
img	input image
return	output image
*/
Mat Dip1::doSomethingThatMyTutorIsGonnaLike(Mat& img){
  
	// TO DO !!!
	Mat inputImage = img;
	Mat outputImage = img;
	
	pixelate(inputImage, outputImage, 8);

	//{yellow, magenta, cyan}
	double color[] = {0.0, -255.0, 0.0};
	stainPicture(outputImage, color);

	return outputImage;
	
}

/* *****************************
  GIVEN FUNCTIONS
***************************** */

// function loads input image, calls processing function, and saves result
/*
fname	path to input image
*/
void Dip1::run(string fname){

	// window names
	string win1 = string ("Original image");
	string win2 = string ("Result");
  
	// some images
	Mat inputImage, outputImage;
  
	// load image
	cout << "load image" << endl;
	inputImage = imread( fname );
	cout << "done" << endl;
	
	// check if image can be loaded
	if (!inputImage.data){
	    cout << "ERROR: Cannot read file " << fname << endl;
	    cout << "Press enter to continue..." << endl;
	    cin.get();
	    exit(-1);
	}

	// show input image
	namedWindow( win1.c_str() );
	imshow( win1.c_str(), inputImage );
	
	// do something (reasonable!)
	outputImage = doSomethingThatMyTutorIsGonnaLike( inputImage );
	
	// show result
	namedWindow( win2.c_str() );
	imshow( win2.c_str(), outputImage );
	
	// save result
	imwrite("result.jpg", outputImage);
	
	// wait a bit
	waitKey(0);

}

// function loads input image and calls the processing functions
// output is tested on "correctness" 
/*
fname	path to input image
*/
void Dip1::test(string fname){

	// some image variables
	Mat inputImage, outputImage;
  
	// load image
	inputImage = imread( fname );

	// check if image can be loaded
	if (!inputImage.data){
	    cout << "ERROR: Cannot read file " << fname << endl;
	    cout << "Continue with pressing enter..." << endl;
	    cin.get();
	    exit(-1);
	}

	// create output
	outputImage = doSomethingThatMyTutorIsGonnaLike( inputImage );
	// test output
	test_doSomethingThatMyTutorIsGonnaLike(inputImage, outputImage);

}

// function loads input image and calls processing function
// output is tested on "correctness" 
/*
inputImage	input image as used by doSomethingThatMyTutorIsGonnaLike()
outputImage	output image as created by doSomethingThatMyTutorIsGonnaLike()
*/
void Dip1::test_doSomethingThatMyTutorIsGonnaLike(Mat& inputImage, Mat& outputImage){

	// ensure that input and output have equal number of channels
	if ( (inputImage.channels() == 3) && (outputImage.channels() == 1) ) //changed and with &&
		cvtColor(inputImage, inputImage, CV_BGR2GRAY);

	// split (multi-channel) image into planes
	vector<Mat> inputPlanes, outputPlanes;
	split( inputImage, inputPlanes );
	split( outputImage, outputPlanes );

	// number of planes (1=grayscale, 3=color)
	int numOfPlanes = inputPlanes.size();

	// calculate and compare image histograms for each plane
	Mat inputHist, outputHist;
	// number of bins
	int histSize = 100;
	float range[] = { 0, 256 } ;
	const float* histRange = { range };
	bool uniform = true; bool accumulate = false;
	double sim = 0;
	for(int p = 0; p < numOfPlanes; p++){
		// calculate histogram
		calcHist( &inputPlanes[p], 1, 0, Mat(), inputHist, 1, &histSize, &histRange, uniform, accumulate );
		calcHist( &outputPlanes[p], 1, 0, Mat(), outputHist, 1, &histSize, &histRange, uniform, accumulate );
		// normalize
		inputHist = inputHist / sum(inputHist).val[0];
		outputHist = outputHist / sum(outputHist).val[0];
		// similarity as histogram intersection
		sim += compareHist(inputHist, outputHist, CV_COMP_INTERSECT);
	}
	sim /= numOfPlanes;

	// check whether images are to similar after transformation
	if (sim >= 0.8)
		cout << "Warning: The input and output image seem to be quite similar (similarity = " << sim << " ). Are you sure your tutor is gonna like your work?" << endl;

}