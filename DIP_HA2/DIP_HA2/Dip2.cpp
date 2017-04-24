//============================================================================
// Name        : Dip2.cpp
// Author      : Ronny Haensch
// Version     : 2.0
// Copyright   : -
// Description : 
//============================================================================

#include "Dip2.h"
#include "list"
#include <cmath>

// convolution in spatial domain
/*
src:     input image
kernel:  filter kernel
return:  convolution result
*/
Mat Dip2::spatialConvolution(Mat& src, Mat& kernel){
   // TO DO !!
	//1. flip kernel
	Mat flipKernel = Mat(kernel.cols, kernel.rows, CV_32FC1);
	for (int x = 0; x < (kernel.cols); x++) {
		for (int y = 0; y < (kernel.rows); y++) {
			kernel.col(kernel.cols - 1 - x).row(kernel.rows - 1 - y).copyTo(flipKernel.col(x).row(y));
		}
	}
	//2. re-center (not needed now)
	//3. multiply and integrate
	Mat outputImage = src.clone();
	Mat mat1 = Mat(3, 3, CV_32F);
	Mat mat2 = Mat(3, 3, CV_32F);
	Rect rect(0, 0, flipKernel.cols, flipKernel.rows);
	for (int x = 0; x < src.cols; x++) {
		rect.x = x - ((kernel.cols - 1)/2);
		for (int y = 0; y < src.rows; y++) {
			rect.y = y - ((kernel.rows - 1)/2);
			//if border -> new pixel = original pixel
			//else -> new pixel = result convolution
			if (x < (kernel.cols - 1) / 2 || x >= (src.cols - (kernel.cols - 1) / 2) ||
				y < (kernel.rows - 1) / 2 || y >= (src.rows - (kernel.rows - 1) / 2)) {
				src.col(x).row(y).copyTo(outputImage.col(x).row(y));
			} else {
				mat1 = src(rect);
				mat2 = mat1.mul(flipKernel);
				outputImage.col(x).row(y) = mean(mat2)*(kernel.cols*kernel.rows);
			};
		}
		cout << ".";
	}
	cout << endl;
   return outputImage;
}

// build spatial kernel
/*
src:     input image
kernel:  filter kernel
return:  convolution result
*/
Mat buildSpatialKernel(int kSize) {

	Mat kernel = Mat::ones(kSize, kSize, CV_32F);
	int anz = 0;
	int powX = 0;
	int powY = 0;
	
	for (int x = 0; x < kSize; x++) {
		if (x >((kSize - 1) / 2)) {
			powX = x - ((kSize - 1) / 2) - 1;
		} else {
			powX = x;
		}
		for (int y = 0; y < kSize; y++) {
			if (y >((kSize - 1) / 2)) {
				powY = y - ((kSize - 1) / 2) - 1;
			} else {
				powY = y;
			}
			//kernel.col(x).row(y) = kernel.col(x).row(y) / (kSize*kSize);
			kernel.col(x).row(y) = kernel.at<float>(Point(x, y)) * pow(2, powX) * pow(2, powY);
			anz += pow(2, powX) * pow(2, powY);
		}
	}
	for (int x = 0; x < kSize; x++) {
		for (int y = 0; y < kSize; y++) {
			//kernel.col(x).row(y) = kernel.col(x).row(y) / (kSize*kSize);
			kernel.col(x).row(y) = kernel.col(x).row(y) / anz;
		}
	}
	return kernel;
}

// build radiometric kernel
/*
src:     input image
anchor:  anchor for kernel
kSize:	 size for each side of the kernel
sigma:	 sigma value for the kernel
return:  radiometric kernel
*/
Mat buildRadiometricKernel(Mat& src, Point2d anchor, int kSize, double sigma) {
	//proof for bound collision
	if (anchor.x - ((kSize - 1) / 2) < 0 ||
		anchor.y - ((kSize - 1) / 2) < 0 ||
		anchor.x + ((kSize - 1) / 2) +1 >= src.cols ||
		anchor.y + ((kSize - 1) / 2) +1 >= src.rows) {
		return Mat::ones(kSize, kSize, CV_32F);
	}
	//get area from src
	Rect rect(anchor.x - ((kSize - 1) / 2), anchor.y - ((kSize - 1) / 2), kSize, kSize);
	Mat srcArea = src(rect);
	Mat rKernel = Mat::ones(kSize, kSize, CV_32F);
	float tmp = 0;
	const float pi = 3.14159265358979323846;
	float rad = (1 / (2 * pi*(sigma*sigma)));
	float pq = 0;
	float srcWeight = 0;
	float anchorWeight = src.at<float>(Point((kSize - 1) / 2, (kSize - 1) / 2));
	//process the formula for radiometric weight
	for (int x = 0; x < kSize; x++) {
		for (int y = 0; y < kSize; y++) {
			srcWeight = srcArea.at<float>(Point(x, y));
			pq = (anchorWeight - srcWeight);
			rKernel.col(x).row(y) = rad * exp(-((pq*pq)/(2*(sigma*sigma))));
			tmp += rKernel.at<float>(Point(x,y));
		}
	}
	//normalize
	tmp = 1 / tmp;
	for (int x = 0; x < kSize; x++) {
		for (int y = 0; y < kSize; y++) {
			rKernel.col(x).row(y) = rKernel.at<float>(Point(x, y)) * tmp;
		}
	}
	return rKernel;
}

/*
src:     input image
bKernel: bilateral kernel
anchor:  anchor for kernel
return:  bilateral weight
*/
float getBilateralWeight(Mat& src, Mat& bKernel, Point2i anchor) {
	//proof for bound collision
	if (anchor.x - ((bKernel.cols - 1) / 2) < 0 ||
		anchor.y - ((bKernel.rows - 1) / 2) < 0 ||
		anchor.x + ((bKernel.cols - 1) / 2) + 1 >= src.cols ||
		anchor.y + ((bKernel.rows - 1) / 2) + 1 >= src.rows) {
		return src.at<float>(anchor);
	}
	//get area from src
	Rect rect(anchor.x - ((bKernel.cols - 1) / 2), anchor.y - ((bKernel.rows - 1) / 2), bKernel.cols, bKernel.rows);
	Mat srcArea = src(rect);
	//multiply filter to area in src
	Mat tmp = srcArea.mul(bKernel);
	//add individual weights to bilateral weight
	float bilateralWeight = sum(tmp).val[0];
	return bilateralWeight;
}

// get median of a Matrix
/*
src:     input image
anchor:  anchor for kernel
kSize:	 size for each side of the kernel
return:  median
*/
float getMedian(Mat& src, Point2i anchor, int size) {

	Point2i medianPixel = anchor;
	Mat kernel = Mat::zeros(size, size, CV_8UC2);
	list<float> list;
	float* array = new float[size*size];
	int iterator = 0;
	
	for (int x = 0; x < size; x++) {
		if (anchor.x - ((size - 1) / 2) + x < 0) {
			return src.at<float>(anchor);
			break;
		} else if (anchor.x - ((size - 1) / 2) + x >= src.cols) {
			return src.at<float>(anchor);
			break;
		} else {
			for (int y = 0; y < size; y++) {
				if (anchor.y - ((size - 1) / 2) + y < 0) {
					return src.at<float>(anchor);
					break;
				} else if (anchor.y - ((size - 1) / 2) + y >= src.rows) {
					return src.at<float>(anchor);
					break;
				} else {
					array[iterator] = src.at<float>(Point(anchor.x - ((size - 1) / 2) + x, anchor.y - ((size - 1) / 2) + y));
					iterator++;
				}
			}
		}
	}

	//int elements = sizeof(array) / sizeof(array[0]);
	sort(array, array + iterator);
	return array[(iterator - 1) / 2];
}

// the average filter
// HINT: you might want to use Dip2::spatialConvolution(...) within this function
/*
src:     input image
kSize:   window size used by local average
return:  filtered image
*/
Mat Dip2::averageFilter(Mat& src, int kSize){
   // TO DO !!
	if (kSize % 2 == 0) {
		cout << "bad filter-size. pls choose odd-number!" << endl;
		return src.clone();
	}
	Mat kernel = Mat::ones(kSize, kSize, CV_32F);
	for (int x = 0; x < kSize; x++) {
		for (int y = 0; y < kSize; y++) {
			kernel.col(x).row(y) = kernel.col(x).row(y) / (kSize*kSize);
		}
	}
	
	Mat outputImage = spatialConvolution(src, kernel);
   return outputImage;
}

// the median filter
/*
src:     input image
kSize:   window size used by median operation
return:  filtered image
*/
Mat Dip2::medianFilter(Mat& src, int kSize){
   // TO DO !!
	if (kSize % 2 == 0) {
		cout << "bad filter-size. pls choose odd-number!" << endl;
		return src.clone();
	}
	Mat outputImage = src.clone();
	Point2i anchor;
	for (int x = 0; x < src.cols; x++) {
		anchor.x = x;
		for (int y = 0; y < src.rows; y++) {
			anchor.y = y;
			outputImage.at<float>(Point(x, y)) = getMedian(src, anchor, kSize);
		}
		cout << ".";
	}
	cout << endl;
	return outputImage;
}

// the bilateral filter
/*
src:     input image
kSize:   size of the kernel --> used to compute std-dev of spatial kernel
sigma:   standard-deviation of the radiometric kernel
return:  filtered image
*/
Mat Dip2::bilateralFilter(Mat& src, int kSize, double sigma){
	if (kSize % 2 == 0) {
		cout << "bad filter-size. pls choose odd-number!" << endl;
		return src.clone();
	}

	Mat outputImage = src.clone();
	//1. build spatial kernel according to kSize
	Mat sKernel = buildSpatialKernel(kSize);
	Mat flipKernel = Mat(kSize, kSize, CV_32FC1);
	//for each pixel in src
	for (int i = 0; i < src.cols; i++) {
		for (int j = 0; j < src.rows; j++) {

			//2. build radiometric kernel according to sigma
			Mat rKernel = buildRadiometricKernel(src, Point(i,j), kSize, sigma);
			//3. build bilateral kernel from spatial and radiometric kernel (multiply and normalize)
			//bilateral kernel not normalized
			Mat bKernel = sKernel.mul(rKernel);
			//get factor for normalization
			float tmp = sum(bKernel).val[0];
			tmp = 1 / tmp;
			//normalize
			for (int x = 0; x < kSize; x++) {
				for (int y = 0; y < kSize; y++) {
					bKernel.col(x).row(y) = bKernel.col(x).row(y) * tmp;
				}
			}
			//4. flip kernel
			for (int x = 0; x < (bKernel.cols); x++) {
				for (int y = 0; y < (bKernel.rows); y++) {
					bKernel.col(bKernel.cols - 1 - x).row(bKernel.rows - 1 - y).copyTo(flipKernel.col(x).row(y));
				}
			}
			bKernel = flipKernel;
			//5. execute filter to src
			outputImage.col(i).row(j) = getBilateralWeight(src, bKernel, Point(i,j));
		}
		cout << ".";
	}
	cout << endl;
	//cv::bilateralFilter(src, outputImage, kSize, sigma, sigma, BORDER_CONSTANT);
    return outputImage;

}

// the non-local means filter
/*
src:   		input image
searchSize: size of search region
sigma: 		Optional parameter for weighting function
return:  	filtered image
*/
Mat Dip2::nlmFilter(Mat& src, int searchSize, double sigma){
  
    return src.clone();

}

/* *****************************
  GIVEN FUNCTIONS
***************************** */

// function loads input image, calls processing function, and saves result
void Dip2::run(void){

   // load images as grayscale
	cout << "load images" << endl;
	Mat noise1 = imread("noiseType_1.jpg", 0);
   if (!noise1.data){
	   cerr << "noiseType_1.jpg not found" << endl;
      cout << "Press enter to exit"  << endl;
      cin.get();
	   exit(-3);
	}
   noise1.convertTo(noise1, CV_32FC1);
	Mat noise2 = imread("noiseType_2.jpg",0);
	if (!noise2.data){
	   cerr << "noiseType_2.jpg not found" << endl;
      cout << "Press enter to exit"  << endl;
      cin.get();
	   exit(-3);
	}
   noise2.convertTo(noise2, CV_32FC1);
	cout << "done" << endl;
	  
   // apply noise reduction
	// TO DO !!!
	// ==> Choose appropriate noise reduction technique with appropriate parameters
	// ==> "average" or "median"? Why?
	// ==> try also "bilateral" (and if implemented "nlm")
	cout << "reduce noise" << endl;
	Mat restorated1 = noiseReduction(noise1, "median", 3, 50);
	Mat restorated2 = noiseReduction(noise2, "bilateral", 3, 50);
	cout << "done" << endl;
	  
	// save images
	cout << "save results" << endl;
	imwrite("restorated1.jpg", restorated1);
	imwrite("restorated2.jpg", restorated2);
	cout << "done" << endl;

}

// noise reduction
/*
src:     input image
method:  name of noise reduction method that shall be performed
	     "average" ==> moving average
         "median" ==> median filter
         "bilateral" ==> bilateral filter
         "nlm" ==> non-local means filter
kSize:   (spatial) kernel size
param:   if method == "bilateral", standard-deviation of radiometric kernel; if method == "nlm", (optional) parameter for similarity function
         can be ignored otherwise (default value = 0)
return:  output image
*/
Mat Dip2::noiseReduction(Mat& src, string method, int kSize, double param){

   // apply moving average filter
   if (method.compare("average") == 0){
      return averageFilter(src, kSize);
   }
   // apply median filter
   if (method.compare("median") == 0){
      return medianFilter(src, kSize);
   }
   // apply bilateral filter
   if (method.compare("bilateral") == 0){
      return bilateralFilter(src, kSize, param);
   }
   // apply adaptive average filter
   if (method.compare("nlm") == 0){
      return nlmFilter(src, kSize, param);
   }

   // if none of above, throw warning and return copy of original
   cout << "WARNING: Unknown filtering method! Returning original" << endl;
   cout << "Press enter to continue"  << endl;
   cin.get();
   return src.clone();

}

// generates and saves different noisy versions of input image
/*
fname:   path to the input image
*/
void Dip2::generateNoisyImages(string fname){
 
   // load image, force gray-scale
   cout << "load original image" << endl;
   Mat img = imread(fname, 0);
   if (!img.data){
      cerr << "ERROR: file " << fname << " not found" << endl;
      cout << "Press enter to exit"  << endl;
      cin.get();
      exit(-3);
   }

   // convert to floating point precision
   img.convertTo(img,CV_32FC1);
   cout << "done" << endl;

   // save original
   imwrite("original.jpg", img);
	  
   // generate images with different types of noise
   cout << "generate noisy images" << endl;

   // some temporary images
   Mat tmp1(img.rows, img.cols, CV_32FC1);
   Mat tmp2(img.rows, img.cols, CV_32FC1);
   // first noise operation
   float noiseLevel = 0.15;
   randu(tmp1, 0, 1);
   threshold(tmp1, tmp2, noiseLevel, 1, CV_THRESH_BINARY);
   multiply(tmp2,img,tmp2);
   threshold(tmp1, tmp1, 1-noiseLevel, 1, CV_THRESH_BINARY);
   tmp1 *= 255;
   tmp1 = tmp2 + tmp1;
   threshold(tmp1, tmp1, 255, 255, CV_THRESH_TRUNC);
   // save image
   imwrite("noiseType_1.jpg", tmp1);
    
   // second noise operation
   noiseLevel = 50;
   randn(tmp1, 0, noiseLevel);
   tmp1 = img + tmp1;
   threshold(tmp1,tmp1,255,255,CV_THRESH_TRUNC);
   threshold(tmp1,tmp1,0,0,CV_THRESH_TOZERO);
   // save image
   imwrite("noiseType_2.jpg", tmp1);

	cout << "done" << endl;
	cout << "Please run now: dip2 restorate" << endl;

}

// function calls some basic testing routines to test individual functions for correctness
void Dip2::test(void){

	test_spatialConvolution();
   test_averageFilter();
   test_medianFilter();

   cout << "Press enter to continue"  << endl;
   cin.get();

}

// checks basic properties of the convolution result
void Dip2::test_spatialConvolution(void){

   Mat input = Mat::ones(9,9, CV_32FC1);
   input.at<float>(4,4) = 255;
   Mat kernel = Mat(3,3, CV_32FC1, 1./9.);

   Mat output = spatialConvolution(input, kernel);
   
   if ( (input.cols != output.cols) || (input.rows != output.rows) ){
      cout << "ERROR: Dip2::spatialConvolution(): input.size != output.size --> Wrong border handling?" << endl;
      return;
   }
  if ( (sum(output.row(0) < 0).val[0] > 0) ||
           (sum(output.row(0) > 255).val[0] > 0) ||
           (sum(output.row(8) < 0).val[0] > 0) ||
           (sum(output.row(8) > 255).val[0] > 0) ||
           (sum(output.col(0) < 0).val[0] > 0) ||
           (sum(output.col(0) > 255).val[0] > 0) ||
           (sum(output.col(8) < 0).val[0] > 0) ||
           (sum(output.col(8) > 255).val[0] > 0) ){
         cout << "ERROR: Dip2::spatialConvolution(): Border of convolution result contains too large/small values --> Wrong border handling?" << endl;
         return;
   }else{
      if ( (sum(output < 0).val[0] > 0) ||
         (sum(output > 255).val[0] > 0) ){
            cout << "ERROR: Dip2::spatialConvolution(): Convolution result contains too large/small values!" << endl;
            return;
      }
   }
   float ref[9][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                      {0, 1, 1, 1, 1, 1, 1, 1, 0},
                      {0, 1, 1, 1, 1, 1, 1, 1, 0},
                      {0, 1, 1, (8+255)/9., (8+255)/9., (8+255)/9., 1, 1, 0},
                      {0, 1, 1, (8+255)/9., (8+255)/9., (8+255)/9., 1, 1, 0},
                      {0, 1, 1, (8+255)/9., (8+255)/9., (8+255)/9., 1, 1, 0},
                      {0, 1, 1, 1, 1, 1, 1, 1, 0},
                      {0, 1, 1, 1, 1, 1, 1, 1, 0},
                      {0, 0, 0, 0, 0, 0, 0, 0, 0}};
   for(int y=1; y<8; y++){
      for(int x=1; x<8; x++){
         if (abs(output.at<float>(y,x) - ref[y][x]) > 0.0001){
            cout << "ERROR: Dip2::spatialConvolution(): Convolution result contains wrong values!" << endl;
            return;
         }
      }
   }
   input.setTo(0);
   input.at<float>(4,4) = 255;
   kernel.setTo(0);
   kernel.at<float>(0,0) = -1;
   output = spatialConvolution(input, kernel);
   if ( abs(output.at<float>(5,5) + 255.) < 0.0001 ){
      cout << "ERROR: Dip2::spatialConvolution(): Is filter kernel \"flipped\" during convolution? (Check lecture/exercise slides)" << endl;
      return;
   }
   if ( ( abs(output.at<float>(2,2) + 255.) < 0.0001 ) || ( abs(output.at<float>(4,4) + 255.) < 0.0001 ) ){
      cout << "ERROR: Dip2::spatialConvolution(): Is anchor point of convolution the centre of the filter kernel? (Check lecture/exercise slides)" << endl;
      return;
   }
   cout << "Message: Dip2::spatialConvolution() seems to be correct" << endl;
}

// checks basic properties of the filtering result
void Dip2::test_averageFilter(void){

   Mat input = Mat::ones(9,9, CV_32FC1);
   input.at<float>(4,4) = 255;

   Mat output = averageFilter(input, 3);
   
   if ( (input.cols != output.cols) || (input.rows != output.rows) ){
      cout << "ERROR: Dip2::averageFilter(): input.size != output.size --> Wrong border handling?" << endl;
      return;
   }
  if ( (sum(output.row(0) < 0).val[0] > 0) ||
           (sum(output.row(0) > 255).val[0] > 0) ||
           (sum(output.row(8) < 0).val[0] > 0) ||
           (sum(output.row(8) > 255).val[0] > 0) ||
           (sum(output.col(0) < 0).val[0] > 0) ||
           (sum(output.col(0) > 255).val[0] > 0) ||
           (sum(output.col(8) < 0).val[0] > 0) ||
           (sum(output.col(8) > 255).val[0] > 0) ){
         cout << "ERROR: Dip2::averageFilter(): Border of result contains too large/small values --> Wrong border handling?" << endl;
         return;
   }else{
      if ( (sum(output < 0).val[0] > 0) ||
         (sum(output > 255).val[0] > 0) ){
            cout << "ERROR: Dip2::averageFilter(): Result contains too large/small values!" << endl;
            return;
      }
   }
   float ref[9][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                      {0, 1, 1, 1, 1, 1, 1, 1, 0},
                      {0, 1, 1, 1, 1, 1, 1, 1, 0},
                      {0, 1, 1, (8+255)/9., (8+255)/9., (8+255)/9., 1, 1, 0},
                      {0, 1, 1, (8+255)/9., (8+255)/9., (8+255)/9., 1, 1, 0},
                      {0, 1, 1, (8+255)/9., (8+255)/9., (8+255)/9., 1, 1, 0},
                      {0, 1, 1, 1, 1, 1, 1, 1, 0},
                      {0, 1, 1, 1, 1, 1, 1, 1, 0},
                      {0, 0, 0, 0, 0, 0, 0, 0, 0}};
   for(int y=1; y<8; y++){
      for(int x=1; x<8; x++){
         if (abs(output.at<float>(y,x) - ref[y][x]) > 0.0001){
            cout << "ERROR: Dip2::averageFilter(): Result contains wrong values!" << endl;
            return;
         }
      }
   }
   cout << "Message: Dip2::averageFilter() seems to be correct" << endl;
}

// checks basic properties of the filtering result
void Dip2::test_medianFilter(void){

   Mat input = Mat::ones(9,9, CV_32FC1);
   input.at<float>(4,4) = 255;

   Mat output = medianFilter(input, 3);
   
   if ( (input.cols != output.cols) || (input.rows != output.rows) ){
      cout << "ERROR: Dip2::medianFilter(): input.size != output.size --> Wrong border handling?" << endl;
      return;
   }
  if ( (sum(output.row(0) < 0).val[0] > 0) ||
           (sum(output.row(0) > 255).val[0] > 0) ||
           (sum(output.row(8) < 0).val[0] > 0) ||
           (sum(output.row(8) > 255).val[0] > 0) ||
           (sum(output.col(0) < 0).val[0] > 0) ||
           (sum(output.col(0) > 255).val[0] > 0) ||
           (sum(output.col(8) < 0).val[0] > 0) ||
           (sum(output.col(8) > 255).val[0] > 0) ){
         cout << "ERROR: Dip2::medianFilter(): Border of result contains too large/small values --> Wrong border handling?" << endl;
         return;
   }else{
      if ( (sum(output < 0).val[0] > 0) ||
         (sum(output > 255).val[0] > 0) ){
            cout << "ERROR: Dip2::medianFilter(): Result contains too large/small values!" << endl;
            return;
      }
   }
   for(int y=1; y<8; y++){
      for(int x=1; x<8; x++){
         if (abs(output.at<float>(y,x) - 1.) > 0.0001){
            cout << "ERROR: Dip2::medianFilter(): Result contains wrong values!" << endl;
            return;
         }
      }
   }
   cout << "Message: Dip2::medianFilter() seems to be correct" << endl;

}
