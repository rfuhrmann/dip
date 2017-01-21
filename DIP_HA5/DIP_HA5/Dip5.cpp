//============================================================================
// Name        : Dip5.cpp
// Author      : Ronny Haensch
// Version     : 2.0
// Copyright   : -
// Description : 
//============================================================================

#define _USE_MATH_DEFINES
#include "Dip5.h"
#include <math.h>

// uses structure tensor to define interest points (foerstner)
void Dip5::getInterestPoints(Mat& img, double sigma, vector<KeyPoint>& points){
	// TO DO !!!
	Mat fstDevKernelX = createFstDevKernel(sigma);
	Mat fstDevKernelY = fstDevKernelX.t();

	//1. gradients in x and y direction
	Mat gradX, gradY, gradXY;
	filter2D(img, gradX, -1, fstDevKernelX);// , Point(-1, -1), NULL, 4);
	filter2D(img, gradY, -1, fstDevKernelY);// , Point(-1, -1), NULL, 4);

	//2. gx*gx, gy*gy, gx*gy
	multiply(gradX, gradY, gradXY);
	pow(gradX, 2, gradX);
	pow(gradY, 2, gradY);
	

	//3. average (gaussian window)
	GaussianBlur(gradXY, gradXY, Size(3, 3), sigma);
	//showImage(gradXY, "Gradients: gx * gy", 0, true, false);
	
	//4. trace of structure tensor
	Mat trace, isotropy;
	add(gradX, gradY, trace);
	//showImage(trace, "Trace", 0, true, false);
	
	//5. determinant of structure tensor
	Mat determinant;
	multiply(gradX, gradY, determinant);
	//showImage(determinant, "Determinant", 0, true, false);

	//6. weight
	Mat weight;
	divide(determinant, trace, weight);
	////showImage(weight, "Weight", 0, true, false);

	//7. weight non-max suppression
	float avgWeight = 0.5*mean(weight)[0];
	weight = nonMaxSuppression(weight);
	//showImage(weight, "Weight - nonMaxSuppression", 0, true, false);

	//8. weight thresholding
	threshold(weight, weight, avgWeight, 1, THRESH_TOZERO);
	//showImage(weight, "Weight - threshold", 0, true, false);

	//9. isotropy
	pow(trace, 2, trace);
	divide(4 * determinant, trace, isotropy);
	//showImage(isotropy, "Isotropy", 0, true, false);

	//10. isotropy non-max suppression
	isotropy = nonMaxSuppression(isotropy);
	//showImage(isotropy, "Isotropy - nonMaxSuppression", 0, true, false);

	//11. isotropy thresholding
	//12. kepoints found
	const float minWeight = 0.01;
	const float minIsotropy = 0.5;
	for (int x = 0; x < img.rows; ++x) {
		for (int y = 0; y < img.cols; ++y) {
			if ((isotropy.at<float>(x, y) >= minIsotropy)
				&& (weight.at<float>(x, y) >= minWeight))
			{
				KeyPoint kp;
				kp.pt.x = y;
				kp.pt.y = x;
				kp.size = 3;
				kp.response = weight.at<float>(x, y);
				cout << "kp" << kp.response<< endl;
				points.push_back(kp);
			}
		}
	}
}

// creates kernel representing fst derivative of a Gaussian kernel in x-direction
/*
sigma	standard deviation of the Gaussian kernel
return	the calculated kernel
*/
Mat Dip5::createFstDevKernel(double sigma){
	// TO DO !!!
	sigma = this->sigma;
	int kSize = (int)ceil(3 * this->sigma);
	kSize += 1 - kSize % 2;
	//int kSize = (((int)ceil(sigma * 3.0)) * 2) + 1;
	Mat gaussianKernel = getGaussianKernel(kSize, sigma, CV_32FC1);
	gaussianKernel = gaussianKernel * gaussianKernel.t();
	Mat fstKernel = Mat::ones(kSize, kSize, CV_32FC1);
	for (int x = 0; x<kSize; x++) {
		for (int y = 0; y<kSize; y++) {
			/*int rx = x - kSize / 2;*/
			fstKernel.at<float>(x, y) = -(x - kSize / 2 )*gaussianKernel.at<float>(x, y) / (sigma * sigma);
		}
	}
	//Sobel
	//fstKernel.at<float>(0, 0) = 1;
	//fstKernel.at<float>(1, 0) = 2;
	//fstKernel.at<float>(2, 0) = 1;
	//fstKernel.at<float>(0, 1) = 0;
	//fstKernel.at<float>(1, 1) = 0;
	//fstKernel.at<float>(2, 1) = 0;
	//fstKernel.at<float>(0, 2) = -1;
	//fstKernel.at<float>(1, 2) = -2;
	//fstKernel.at<float>(2, 2) = -1;

	//cout << fstKernel << endl;
	return fstKernel;
}

/* *****************************
  GIVEN FUNCTIONS
***************************** */

// function calls processing function
/*
in		:  input image
points	:	detected keypoints
*/
void Dip5::run(Mat& in, vector<KeyPoint>& points){
   this->getInterestPoints(in, this->sigma, points);
}

// non-maxima suppression
// if any of the pixel at the 4-neighborhood is greater than current pixel, set it to zero
Mat Dip5::nonMaxSuppression(Mat& img){

	Mat out = img.clone();
	
	for(int x=1; x<out.cols-1; x++){
		for(int y=1; y<out.rows-1; y++){
			if ( img.at<float>(y-1, x) >= img.at<float>(y, x) ){
				out.at<float>(y, x) = 0;
				continue;
			}
			if ( img.at<float>(y, x-1) >= img.at<float>(y, x) ){
				out.at<float>(y, x) = 0;
				continue;
			}
			if ( img.at<float>(y, x+1) >= img.at<float>(y, x) ){
				out.at<float>(y, x) = 0;
				continue;
			}
			if ( img.at<float>( y+1, x) >= img.at<float>(y, x) ){
				out.at<float>(y, x) = 0;
				continue;
			}
		}
	}
	return out;
}

// Function displays image (after proper normalization)
/*
win   :  Window name
img   :  Image that shall be displayed
cut   :  whether to cut or scale values outside of [0,255] range
*/
void Dip5::showImage(Mat& img, const char* win, int wait, bool show, bool save){
  
    Mat aux = img.clone();

    // scale and convert
    if (img.channels() == 1)
		normalize(aux, aux, 0, 255, CV_MINMAX);
		aux.convertTo(aux, CV_8UC1);
    // show
    if (show){
      imshow( win, aux);
      waitKey(wait);
    }
    // save
    if (save)
      imwrite( (string(win)+string(".png")).c_str(), aux);
}
