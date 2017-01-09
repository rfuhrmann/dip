//============================================================================
// Name        : Dip4.cpp
// Author      : Ronny Haensch
// Version     : 2.0
// Copyright   : -
// Description : 
//============================================================================

#include "Dip4.h"

// Performes a circular shift in (dx,dy) direction
/*
in       :  input matrix
dx       :  shift in x-direction
dy       :  shift in y-direction
return   :  circular shifted matrix
*/
Mat Dip4::circShift(Mat& in, int dx, int dy){
	// TO DO !!!
	Mat out = in.clone();
	int currX = 0;
	int currY = 0;
	for (int x = 0; x < in.cols; x++) {
		currX = x - dx;
		//if walk through left border
		if (currX < 0) {
			currX = in.cols + (currX%in.cols);
		}
		//if walk through right border
		else if (currX >= in.cols) {
			currX = currX % in.cols;
		}
		for (int y = 0; y < in.rows; y++) {
			currY = y - dy;
			//if walk through top border
			if (currY < 0) {
				currY = in.rows + (currY%in.rows);
			}
			//if walk through bottom border
			else if (currY >= in.rows) {
				currY = currY % in.rows;
			}
			in.col(currX).row(currY).copyTo(out.col(x).row(y));
		}
	}
	return out;
}

// Function applies inverse filter to restorate a degraded image
/*
degraded :  degraded input image
filter   :  filter which caused degradation
return   :  restorated output image
*/
Mat Dip4::inverseFilter(Mat& degraded, Mat& filter){
	// TO DO !!!

	//Rect rect(4, 4, 4, 4);
	////cout << "rect1" << fKernel(rect) << endl;
	////for (int y = 0; y < fKernel.rows; y++) {
	////	for (int x = 0; x < fKernel.cols; x++) {
	////		if (abs(fKernel.at<float>(Point(x, y))) < t) {
	////			fKernel.at<float>(Point(x, y)) =  1 / t;
	////		}
	////		else {
	////			fKernel.col(x).row(y) = 1 / (fKernel.at<float>(Point(x, y)));
	////		}
	////	}
	////}
	////cout << "rect2" << fKernel(rect) << endl;

	Mat fImage, fKernel, outputImage;

	//copy kernel to large matrix (the size of input image)
	Rect rect(Point(0, 0), filter.size());
	fKernel = Mat::zeros(degraded.size(), degraded.type());
	filter.copyTo(fKernel(rect));

	//circ shift the kernel 
	fKernel = circShift(fKernel, -floor(filter.cols / 2), -floor(filter.rows / 2));

	//convert to frequency domains
	dft(degraded, fImage, DFT_COMPLEX_OUTPUT);
	dft(fKernel, fKernel, DFT_COMPLEX_OUTPUT);

	Mat fPlanes[2];
	split(fKernel, fPlanes);

	// T = epsilon * maxMagnitude(fKernel)
	const double epsilon = 0.05f;
	Mat Magnitude;
	magnitude(fPlanes[0], fPlanes[1], Magnitude);
	double maxMagnitude;
	minMaxIdx(abs(Magnitude), 0, &maxMagnitude, 0);
	const double threshold = epsilon * maxMagnitude;
	

	//Replace the inverse filter 1/Pi by Qi
	for (int y = 0; y < fKernel.rows; y++) {
		for (int x = 0; x < fKernel.cols; x++) {
			if (norm(fKernel.at<Vec2f>(Point(x, y))) < threshold) {
				fKernel.at<Vec2f>(Point(x, y)) = threshold;
			}
			complex<float> p(fKernel.at<Vec2f>(Point(x,y))[0], fKernel.at<Vec2f>(Point(x, y))[1]);
			p = 1.f / p;
			fKernel.at<Vec2f>(Point(x, y))[0] = p.real();
			fKernel.at<Vec2f>(Point(x, y))[1] = p.imag();
			//divide(1.f, fPlanes[0], fPlanes[0]);
			//divide(1.f, fPlanes[1], fPlanes[1]);
			//fKernel.at<Vec2f>(Point(x,y))[0]
		}
	}

	//multiplication in freq domains
	mulSpectrums(fImage, fKernel, fImage, 0);

	//convert to spatial domain
	dft(fImage, outputImage, DFT_INVERSE | DFT_REAL_OUTPUT + DFT_SCALE);

	return outputImage;

}

// Function applies wiener filter to restorate a degraded image
/*
degraded :  degraded input image
filter   :  filter which caused degradation
snr      :  signal to noise ratio of the input image
return   :   restorated output image
*/
Mat Dip4::wienerFilter(Mat& degraded, Mat& filter, double snr){
	// TO DO !!!
	Mat fImage, fKernel, outputImage;

	//copy kernel to large matrix (the size of input image)
	Rect rect(Point(0, 0), filter.size());
	fKernel = Mat::zeros(degraded.size(), degraded.type());
	filter.copyTo(fKernel(rect));

	//perform circ shift on kernel 
	fKernel = circShift(fKernel, -filter.cols / 2, -filter.rows / 2);

	//convert to frequency domains
	dft(degraded, fImage, DFT_COMPLEX_OUTPUT);
	dft(fKernel, fKernel, DFT_COMPLEX_OUTPUT);

	Mat fPlanes[2];
	split(fKernel, fPlanes);

	//P*
	for (int x = 0; x < fKernel.cols; x++) {
		for (int y = 0; y < fKernel.rows; y++) {
			//imaginary part
			fPlanes[1].at<float>(Point(x,y)) = -1 * fPlanes[1].at<float>(Point(x,y));
		}
	}

	Mat fMagnitude;
	magnitude(fPlanes[0], fPlanes[1], fMagnitude);
	multiply(fMagnitude, fMagnitude, fMagnitude);

	//Q = P* / ( |P|² + 1/SNR² )
	divide(fPlanes[0], fMagnitude + 1 / (snr * snr), fPlanes[0]);
	divide(fPlanes[1], fMagnitude + 1 / (snr * snr), fPlanes[1]);
	merge(fPlanes, 2, fKernel);

	//multiplication in freq domains
	mulSpectrums(fImage, fKernel, fImage, 0);

	//convert to spatial domain
	dft(fImage, outputImage, DFT_INVERSE | DFT_REAL_OUTPUT + DFT_SCALE);

	//threshold(outputImage, outputImage, 255, 0, THRESH_TRUNC);
	threshold(outputImage, outputImage, 255, 255, CV_THRESH_TRUNC);
	threshold(outputImage, outputImage, 0, 0, CV_THRESH_TOZERO);
	return outputImage;
}

/* *****************************
  GIVEN FUNCTIONS
***************************** */

// function calls processing function
/*
in                   :  input image
restorationType     :  integer defining which restoration function is used
kernel               :  kernel used during restoration
snr                  :  signal-to-noise ratio (only used by wieder filter)
return               :  restorated image
*/
Mat Dip4::run(Mat& in, string restorationType, Mat& kernel, double snr){

   if (restorationType.compare("wiener")==0){
      return wienerFilter(in, kernel, snr);
   }else{
      return inverseFilter(in, kernel);
   }

}

// function degrades the given image with gaussian blur and additive gaussian noise
/*
img         :  input image
degradedImg :  degraded output image
filterDev   :  standard deviation of kernel for gaussian blur
snr         :  signal to noise ratio for additive gaussian noise
return      :  the used gaussian kernel
*/
Mat Dip4::degradeImage(Mat& img, Mat& degradedImg, double filterDev, double snr){

    int kSize = round(filterDev*3)*2 - 1;
   
    Mat gaussKernel = getGaussianKernel(kSize, filterDev, CV_32FC1);
    gaussKernel = gaussKernel * gaussKernel.t();

    Mat imgs = img.clone();
    dft( imgs, imgs, CV_DXT_FORWARD, img.rows);
    Mat kernels = Mat::zeros( img.rows, img.cols, CV_32FC1);
    int dx, dy; dx = dy = (kSize-1)/2.;
    for(int i=0; i<kSize; i++) for(int j=0; j<kSize; j++) kernels.at<float>((i - dy + img.rows) % img.rows,(j - dx + img.cols) % img.cols) = gaussKernel.at<float>(i,j);
	dft( kernels, kernels, CV_DXT_FORWARD );
	mulSpectrums( imgs, kernels, imgs, 0 );
	dft( imgs, degradedImg, CV_DXT_INV_SCALE, img.rows );
	
    Mat mean, stddev;
    meanStdDev(img, mean, stddev);

    Mat noise = Mat::zeros(img.rows, img.cols, CV_32FC1);
    randn(noise, 0, stddev.at<double>(0)/snr);
    degradedImg = degradedImg + noise;
    threshold(degradedImg, degradedImg, 255, 255, CV_THRESH_TRUNC);
    threshold(degradedImg, degradedImg, 0, 0, CV_THRESH_TOZERO);

    return gaussKernel;
}

// Function displays image (after proper normalization)
/*
win   :  Window name
img   :  Image that shall be displayed
cut   :  whether to cut or scale values outside of [0,255] range
*/
void Dip4::showImage(const char* win, Mat img, bool cut){

   Mat tmp = img.clone();

   if (tmp.channels() == 1){
      if (cut){
         threshold(tmp, tmp, 255, 255, CV_THRESH_TRUNC);
         threshold(tmp, tmp, 0, 0, CV_THRESH_TOZERO);
      }else
         normalize(tmp, tmp, 0, 255, CV_MINMAX);
         
      tmp.convertTo(tmp, CV_8UC1);
   }else{
      tmp.convertTo(tmp, CV_8UC3);
   }
   imshow(win, tmp);
}

// function calls basic testing routines to test individual functions for correctness
void Dip4::test(void){

   test_circShift();
   cout << "Press enter to continue"  << endl;
   cin.get();

}

void Dip4::test_circShift(void){
   
   Mat in = Mat::zeros(3,3,CV_32FC1);
   in.at<float>(0,0) = 1;
   in.at<float>(0,1) = 2;
   in.at<float>(1,0) = 3;
   in.at<float>(1,1) = 4;
   Mat ref = Mat::zeros(3,3,CV_32FC1);
   ref.at<float>(0,0) = 4;
   ref.at<float>(0,2) = 3;
   ref.at<float>(2,0) = 2;
   ref.at<float>(2,2) = 1;
   
   if (sum((circShift(in, -1, -1) == ref)).val[0]/255 != 9){
      cout << "ERROR: Dip4::circShift(): Result of circshift seems to be wrong!" << endl;
      return;
   }
   cout << "Message: Dip4::circShift() seems to be correct" << endl;
}
