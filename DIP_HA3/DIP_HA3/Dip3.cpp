//============================================================================
// Name    : Dip3.cpp
// Author   : Ronny Haensch
// Version    : 2.0
// Copyright   : -
// Description : 
//============================================================================

#include "Dip3.h"
//#include <Windows.h>

// Generates a gaussian filter kernel of given size
/*
kSize:     kernel size (used to calculate standard deviation)
return:    the generated filter kernel
*/
Mat Dip3::createGaussianKernel(int kSize) {
	float kernelMid = (kSize - 1) / 2;
	float sigma = kSize / 5;

	Mat spatialKernel = Mat::ones(kSize, kSize, CV_32FC1);
	for (int x = 0; x < kSize; x++) {
		for (int y = 0; y < kSize; y++) {
				spatialKernel.at<float>(y, x) = exp(-0.5 * (pow(x - kernelMid, 2) / sigma + pow(y - kernelMid, 2) / sigma));
		}
	}
	spatialKernel = spatialKernel / sum(spatialKernel).val[0];
	return spatialKernel;
}

// Performes a circular shift in (dx,dy) direction
/*
in       input matrix
dx       shift in x-direction
dy       shift in y-direction
return   circular shifted matrix
*/
Mat Dip3::circShift(Mat& in, int dx, int dy) {

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

//Performes a convolution by multiplication in frequency domain
/*
in       input image
kernel   filter kernel
return   output image
*/
Mat Dip3::frequencyConvolution(Mat& in, Mat& kernel) {

	// TO DO !!!
	Mat fKernel = Mat::zeros(in.rows, in.cols, CV_32FC1);
	kernel.copyTo(fKernel.colRange(0, kernel.cols).rowRange(0, kernel.rows));
	//1. circShift Kernel
	fKernel = circShift(fKernel, -((kernel.cols - 1) / 2), -((kernel.rows - 1) / 2));

	//2. build Kernel in size of image
	Mat fImage, outputImage;
	in.copyTo(fImage);

	//3. dft of kernel and image
	dft(fImage, fImage, 0);
	dft(fKernel, fKernel, 0);

	//4. multiply spectrums of kernel and image in frequenzy domain
	mulSpectrums(fImage, fKernel, outputImage, 0);

	//5. inverse dft
	dft(outputImage, outputImage, DFT_INVERSE + DFT_SCALE);

	return outputImage;
}

//Performes a threshold by multiplication in frequency domain
//Set Value of pixel to 0 if lower bound for threshold not reached
/*
image    input image in frequency domain
thresh   magnitude of thresholding
return   fImage in frequenzydomain
*/
Mat threshold(Mat& image, double thresh) {
	Mat fImage = image.clone();
	dft(fImage, fImage, DFT_INVERSE + DFT_SCALE);
	for (int y = 0; y < fImage.rows; y++) {
		for (int x = 0; x < fImage.cols; x++) {
			if (fImage.at<float>(Point(x, y)) <= thresh) { //abs
				fImage.at<float>(Point(x, y)) = 0;
			}
		}
	}
	dft(fImage, fImage, 0);
	return fImage;
}

// Performs UnSharp Masking to enhance fine image structures
/*
in       the input image
type     integer defining how convolution for smoothing operation is done
0 <==> spatial domain; 1 <==> frequency domain; 2 <==> seperable filter; 3 <==> integral image
size     size of used smoothing kernel
thresh   minimal intensity difference to perform operation
scale    scaling of edge enhancement
return   enhanced image
*/
//Mat Dip3::usm(Mat& in, int type, int size, double thresh, double scale) {
//	########## Frequency Domain ##########
//	// some temporary images 
//	Mat tmp(in.rows, in.cols, CV_32FC1);
//
//	// calculate edge enhancement
//
//	// 1: smooth original image
//	//    save result in tmp for subsequent usage
//	switch (type) {
//	case 0:
//		tmp = mySmooth(in, size, 0);
//		break;
//	case 1:
//		tmp = mySmooth(in, size, 1);
//		break;
//	case 2:
//		tmp = mySmooth(in, size, 2);
//		break;
//	case 3:
//		tmp = mySmooth(in, size, 3);
//		break;
//	default:
//		GaussianBlur(in, tmp, Size(floor(size / 2) * 2 + 1, floor(size / 2) * 2 + 1), size / 5., size / 5.);
//	}
//
//	// TO DO !!!
//	/*Smooth : y0 -> y1
//	2. Subtract : y2 = y0 – y1
//	2.5.Scale : ý2 = s * y2
//	3. Add : y3 = y0 + ý2
//	-> Final(green)*/
//	Mat fImage, fIn;
//	//original image
//	dft(in, fIn, 0);
//	//smoothed image
//	fImage = tmp;
//	dft(fImage, fImage, 0);
//
//	//2. subtract
//	subtract(fIn, fImage, fImage);
//
//	//3. threshold
//	//fImage = threshold(fImage, thresh);
//	dft(fImage, fImage, DFT_INVERSE + DFT_SCALE);
//	threshold(fImage, fImage, thresh, 255, CV_THRESH_TOZERO);
//	dft(fImage, fImage, 0);
//	//4. scale
//	fImage = fImage * scale;// (1 / abs(scale*scale));
//
//	//5. add
//	add(fImage, fIn, fImage);
//
//	dft(fImage, fImage, DFT_INVERSE + DFT_SCALE);
//	threshold(fImage, fImage, 255, 255, CV_THRESH_TRUNC);
//	threshold(fImage, fImage, 0, 0, CV_THRESH_TOZERO);
//	return fImage;
//}
Mat Dip3::usm(Mat& in, int type, int size, double thresh, double scale) {

	Mat tmp(in.rows, in.cols, CV_32FC1);

	// calculate edge enhancement
	// 1: smooth original image
	//    save result in tmp for subsequent usage
	switch (type) {
	case 0:
		tmp = mySmooth(in, size, 0);
		break;
	case 1:
		tmp = mySmooth(in, size, 1);
		break;
	case 2:
		tmp = mySmooth(in, size, 2);
		break;
	case 3:
		tmp = mySmooth(in, size, 3);
		break;
	default:
		GaussianBlur(in, tmp, Size(floor(size / 2) * 2 + 1, floor(size / 2) * 2 + 1), size / 5., size / 5.);
	}

	// TO DO !!!
	//smoothed image
	Mat fImage = tmp;
	//dft(fImage, fImage, 0);

	//2. subtract
	subtract(in, tmp, fImage);

	for (int y = 0; y < fImage.rows; y++) {
		for (int x = 0; x < fImage.cols; x++) {
			//3. threshold
			if (abs(fImage.at<float>(Point(x, y)) > thresh)) {
				//4. scale
				fImage.at<float>(Point(x, y)) = fImage.at<float>(Point(x, y))*scale;
			}
		}
	}
	//for (int y = 0; y < fImage.rows; y++) {
	//	for (int x = 0; x < fImage.cols; x++) {
	//		//3. threshold
	//		if (abs(fImage.at<float>(Point(x, y)) <= thresh)) {
	//			fImage.at<float>(Point(x, y)) = 0;
	//		}
	//	}
	//}
	////4. scale
	//fImage = fImage * scale;

	//5. add
	fImage = fImage + in;
	//add(fImage, fIn, fImage);

	return fImage;
}


// Build Matrix when jumping over boarder
// not used, cause result is right, but look bad
/*
src:    input image
kernel:  filter kernel
return:  image matrix for convolution
*/
Mat spatialConvolutionBorderHelper(Mat& src, Mat& kernel, Point2i anchor) {
	Mat out = Mat::zeros(kernel.cols, kernel.rows, CV_32F);
	Vec2b vec;
	for (int y = 0; y < kernel.rows; y++) {
		for (int x = 0; x < kernel.cols; x++) {
			vec = (src.cols + (anchor.x - ((kernel.cols - 1) / 2) + x) % src.cols, src.rows + (anchor.y - ((kernel.rows - 1) / 2) + y) % src.rows);
			src.col(vec[0]).row(vec[1]).copyTo(out.col(x).row(y));
		}
	}
	return out;
}

// convolution in spatial domain
/*
src:    input image
kernel:  filter kernel
return:  convolution result
*/
Mat Dip3::spatialConvolution(Mat& src, Mat& kernel) {

	// Hopefully already DONE
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
	Mat mat1 = Mat(flipKernel.cols, flipKernel.rows, CV_32FC1);
	Mat mat2 = Mat(flipKernel.cols, flipKernel.rows, CV_32FC1);
	Rect rect(0, 0, flipKernel.cols, flipKernel.rows);
	for (int x = 0; x < src.cols; x++) {
		rect.x = x - ((kernel.cols - 1) / 2);
		for (int y = 0; y < src.rows; y++) {
			rect.y = y - ((kernel.rows - 1) / 2);
			//if border -> new pixel = original pixel
			//else -> new pixel = result convolution
			if (x < (kernel.cols - 1) / 2 || x >= (src.cols - (kernel.cols - 1) / 2) ||
				y < (kernel.rows - 1) / 2 || y >= (src.rows - (kernel.rows - 1) / 2)) {
				src.col(x).row(y).copyTo(outputImage.col(x).row(y));
				//for (int a = 0; a < kernel.rows*10; a++) {
				//	for (int b = 0; b < kernel.cols*10; b++) {
				//		float tmp = 123456789 / 12345678987654321;
				//	}
				//}
				//mat1 = spatialConvolutionBorderHelper(src, kernel, Point2i(x,y));
			}
			else {
				mat1 = src(rect);
				mat2 = mat1.mul(flipKernel);
				outputImage.col(x).row(y) = mean(mat2)*(kernel.cols*kernel.rows);
			};
		}
		//cout << ".";
	}
	//cout << endl;
	return outputImage;
}


// convolution in spatial domain by seperable filters
/*
src:    input image
size     size of filter kernel
return:  convolution result
*/
Mat Dip3::seperableFilter(Mat& src, int size) {

	// optional

	return src;

}

// convolution in spatial domain by integral images
/*
src:    input image
size     size of filter kernel
return:  convolution result
*/
Mat Dip3::satFilter(Mat& src, int size) {

	// optional

	return src;

}

/* *****************************
GIVEN FUNCTIONS
***************************** */

// function calls processing function
/*
in       input image
type     integer defining how convolution for smoothing operation is done
0 <==> spatial domain; 1 <==> frequency domain
size     size of used smoothing kernel
thresh   minimal intensity difference to perform operation
scale    scaling of edge enhancement
return   enhanced image
*/
Mat Dip3::run(Mat& in, int smoothType, int size, double thresh, double scale) {

	return usm(in, smoothType, size, thresh, scale);

}


// Performes smoothing operation by convolution
/*
in       input image
size     size of filter kernel
type     how is smoothing performed?
return   smoothed image
*/
Mat Dip3::mySmooth(Mat& in, int size, int type) {

	// create filter kernel
	Mat kernel = createGaussianKernel(size);

	// perform convoltion
	switch (type) {
	case 0: return spatialConvolution(in, kernel);	// 2D spatial convolution
	case 1: return frequencyConvolution(in, kernel);	// 2D convolution via multiplication in frequency domain
	case 2: return seperableFilter(in, size);	// seperable filter
	case 3: return satFilter(in, size);		// integral image
	default: return frequencyConvolution(in, kernel);
	}
}

// function calls basic testing routines to test individual functions for correctness
void Dip3::test(void) {

	test_createGaussianKernel();
	test_circShift();
	test_frequencyConvolution();
	cout << "Press enter to continue" << endl;
	cin.get();

}

void Dip3::test_createGaussianKernel(void) {

	Mat k = createGaussianKernel(11);

	if (abs(sum(k).val[0] - 1) > 0.0001) {
		cout << "ERROR: Dip3::createGaussianKernel(): Sum of all kernel elements is not one!" << endl;
		return;
	}
	if (sum(k >= k.at<float>(5, 5)).val[0] / 255 != 1) {
		cout << "ERROR: Dip3::createGaussianKernel(): Seems like kernel is not centered!" << endl;
		return;
	}
	cout << "Message: Dip3::createGaussianKernel() seems to be correct" << endl;
}

void Dip3::test_circShift(void) {

	Mat in = Mat::zeros(3, 3, CV_32FC1);
	in.at<float>(0, 0) = 1;
	in.at<float>(0, 1) = 2;
	in.at<float>(1, 0) = 3;
	in.at<float>(1, 1) = 4;
	Mat ref = Mat::zeros(3, 3, CV_32FC1);
	ref.at<float>(0, 0) = 4;
	ref.at<float>(0, 2) = 3;
	ref.at<float>(2, 0) = 2;
	ref.at<float>(2, 2) = 1;

	if (sum((circShift(in, -1, -1) == ref)).val[0] / 255 != 9) {
		cout << "ERROR: Dip3::circShift(): Result of circshift seems to be wrong!" << endl;
		return;
	}
	cout << "Message: Dip3::circShift() seems to be correct" << endl;
}

void Dip3::test_frequencyConvolution(void) {

	Mat input = Mat::ones(9, 9, CV_32FC1);
	input.at<float>(4, 4) = 255;
	Mat kernel = Mat(3, 3, CV_32FC1, 1. / 9.);

	Mat output = frequencyConvolution(input, kernel);

	if ((sum(output < 0).val[0] > 0) || (sum(output > 255).val[0] > 0)) {
		cout << "ERROR: Dip3::frequencyConvolution(): Convolution result contains too large/small values!" << endl;
		return;
	}
	float ref[9][9] = { { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 1, 1, 1, 1, 1, 1, 1, 0 },
	{ 0, 1, 1, 1, 1, 1, 1, 1, 0 },
	{ 0, 1, 1, (8 + 255) / 9., (8 + 255) / 9., (8 + 255) / 9., 1, 1, 0 },
	{ 0, 1, 1, (8 + 255) / 9., (8 + 255) / 9., (8 + 255) / 9., 1, 1, 0 },
	{ 0, 1, 1, (8 + 255) / 9., (8 + 255) / 9., (8 + 255) / 9., 1, 1, 0 },
	{ 0, 1, 1, 1, 1, 1, 1, 1, 0 },
	{ 0, 1, 1, 1, 1, 1, 1, 1, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
	for (int y = 1; y<8; y++) {
		for (int x = 1; x<8; x++) {
			if (abs(output.at<float>(y, x) - ref[y][x]) > 0.0001) {
				cout << "ERROR: Dip3::frequencyConvolution(): Convolution result contains wrong values!" << endl;
				return;
			}
		}
	}
	cout << "Message: Dip3::frequencyConvolution() seems to be correct" << endl;
}
