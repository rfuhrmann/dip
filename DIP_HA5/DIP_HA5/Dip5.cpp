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
	cout << "Interest Points calculation" << endl;

	cout << "initialize stuff" << endl;
	Mat gauss = img.clone(); //= createGaussianKernel(2 * sigma);
	GaussianBlur(img, gauss, Size(0,0), sigma, sigma, 4);
	Mat devXK = createFstDevKernel(sigma);
	Mat devYK;
	transpose(devXK, devYK);


	cout << "calculate the x- and y-gradients" << endl;
	//Mat gradX = spatialConvolution(img, devXK);
	//Mat gradY = spatialConvolution(img, devYK);
	Mat gradX, gradY, gradXY;
	filter2D(img, gradX, -1, devXK, Point(-1, -1), NULL, 4);
	filter2D(img, gradY, -1, devYK, Point(-1, -1), NULL, 4);

	multiply(gradX, gradY, gradXY);
	showImage(gradXY, "Gradients: gx * gy", 0, true, false);
	
	pow(gradX, 2, gradX);
	pow(gradY, 2, gradY);
	Mat trace, determinant, weight, isotropy;
	add(gradX, gradY, trace);
	showImage(trace, "Trace", 0, true, false);
	multiply(gradX, gradY, determinant);
	showImage(determinant, "Determinant", 0, true, false);
	divide(determinant, trace, weight);
	showImage(weight, "Weight", 0, true, false);
	nonMaxSuppression(weight);
	showImage(weight, "Weight - nonMaxSuppression", 0, true, false);

	//WEIGHT TRASHOLDING

	multiply(weight, 4, weight);
	pow(trace, 2, trace);
	divide(weight, trace, isotropy);
	showImage(isotropy, "Isotropy", 0, true, false);
	nonMaxSuppression(isotropy);
	showImage(isotropy, "Isotropy - nonMaxSuppression", 0, true, false);

	//ISOTROPY TRASHOLDING

	const float wMin = 70.0; // 0.5 ... 1.5
	const float qMin = 0.5; // 0.5 ... 0.75
	const float wMinActual = wMin * mean(weight)[0];
	for (int x = 0; x < img.rows; ++x) {
		for (int y = 0; y < img.cols; ++y) {
			if ((isotropy.at<float>(x, y) >= qMin)
				&& (weight.at<float>(x, y) >= wMinActual))
			{
				KeyPoint kp;
				kp.pt.x = y;
				kp.pt.y = x;
				kp.angle = 0;
				kp.size = 3;
				kp.response = weight.at<float>(x, y);
				points.push_back(kp);
			}
		}
	}

	//Mat gXX;
	//Mat gYY;
	//Mat gXY;
	//multiply(gradX, gradY, gXY);
	//multiply(gradX, gradX, gXX);
	//multiply(gradY, gradY, gYY);

	////Mat gXXs = spatialConvolution(gXX, gauss);
	////Mat gYYs = spatialConvolution(gYY, gauss);
	////Mat gXYs = spatialConvolution(gXY, gauss);
	//Mat gXXs, gYYs, gXYs;
	//filter2D(gauss, gXXs, -1, gXX, Point(-1, -1), NULL, 4);
	//filter2D(gauss, gYYs, -1, gYY, Point(-1, -1), NULL, 4);
	//filter2D(gauss, gXYs, -1, gXY, Point(-1, -1), NULL, 4);

	//cout << "calculate the trace" << endl;
	//// see DIP05_ST13_interest.pdf page 22
	//Mat trace, traceX, traceY;
	//add(gXXs, gYYs, trace);
	//
	////test
	//pow(gXXs, 2, traceX);
	//pow(gYYs, 2, traceY);
	//add(traceX, traceY, traceX);
	////pow(trace, 2, traceX);
	//showImage(traceX, "Trace", 0, true, false);

	//cout << "calculate the determinant" << endl;
	//// see DIP05_ST13_interest.pdf page 22
	//Mat det;
	//Mat temp;
	//multiply(gXXs, gYYs, det);
	//multiply(gXYs, gXYs, temp);
	//subtract(det, temp, det);

	//showImage(det, "Determinant", 0, true, false);

	//cout << "calculate the weight" << endl;
	//// Strength of gradients in the neighborhood
	//// see DIP05_ST13_interest.pdf page 22
	//Mat w;
	//divide(det, trace, w);

	//showImage(w, "Weight", 0, true, false);

	//Mat wNoMax;
	////nonMaxSuppression(w, wNoMax);
	//nonMaxSuppression(w);
	//w.copyTo(wNoMax);
	//cout << "calculate the isotropy" << endl;
	//// Measures the uniformity of gradient directions in the neighbourhood
	//// see DIP05_ST13_interest.pdf page 22
	//multiply(trace, trace, temp);
	//Mat q;
	//divide(4 * det, temp, q);
	//Mat qNoMax;
	////nonMaxSuppression(q, qNoMax);
	//nonMaxSuppression(temp);
	//temp.copyTo(qNoMax);

	//showImage(qNoMax, "Isotropy", 0, true, false);

	//cout << "do the thresholding" << endl;
	//// XXX why do we have to choose wMin so high?
	//const float wMin = 70.0; // 0.5 ... 1.5
	//const float qMin = 0.5; // 0.5 ... 0.75
	//const float wMinActual = wMin * mean(wNoMax)[0];
	//for (int x = 0; x < img.rows; ++x) {
	//	for (int y = 0; y < img.cols; ++y) {
	//		if ((qNoMax.at<float>(x, y) >= qMin)
	//			&& (wNoMax.at<float>(x, y) >= wMinActual))
	//		{
	//			KeyPoint kp;
	//			kp.pt.x = y;
	//			kp.pt.y = x;
	//			kp.angle = 0;
	//			kp.size = 3;
	//			kp.response = wNoMax.at<float>(x, y);
	//			points.push_back(kp);
	//		}
	//	}
	//}
	//showImage(img, "keypoints", 0, true, false);
}

// creates kernel representing fst derivative of a Gaussian kernel in x-direction
/*
sigma	standard deviation of the Gaussian kernel
return	the calculated kernel
*/
Mat Dip5::createFstDevKernel(double sigma){
	// TO DO !!!
	const int kSize = (((int)ceil(sigma * 3.0)) * 2) + 1;
	Mat kernel(1, kSize, CV_32FC1);

	const int w = kernel.rows;
	const int h = kernel.cols;
	const float sigmaSqr = sigma * sigma;
	const float sigma22 = 2 * sigmaSqr;
	const float sigma4Pi = 2 * M_PI * sigmaSqr * sigmaSqr;
	const float hw = (w / 2);
	const float hh = (h / 2);
	for (int x = 0; x < w; ++x) {
		const float xv = x - hw;
		for (int y = 0; y < h; ++y) {
			const float yv = y - hh;
			const float commonFac = exp(-(xv*xv + yv*yv) / sigma22) / sigma4Pi;

			// see DIP05_ST13_interest.pdf page 10
			if (kernel.type() == CV_32FC2) {
				kernel.at<Vec2f>(x, y)[0] = -xv * commonFac; // partial derivative towards x
				kernel.at<Vec2f>(x, y)[1] = -yv * commonFac; // partial derivative towards y
			}
			else {
				kernel.at<float>(x, y) = -yv * commonFac; // partial derivative towards y
			}
		}
	}
	return kernel;//Mat::zeros(3, 3, CV_32FC1);
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
