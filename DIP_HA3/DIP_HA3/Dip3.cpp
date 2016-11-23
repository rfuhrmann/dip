//============================================================================
// Name    : Dip3.cpp
// Author   : Ronny Haensch
// Version    : 2.0
// Copyright   : -
// Description : 
//============================================================================

#include "Dip3.h"

// Generates a gaussian filter kernel of given size
/*
kSize:     kernel size (used to calculate standard deviation)
return:    the generated filter kernel
*/
Mat Dip3::createGaussianKernel(int kSize) {

	// TO DO !!!
	if (kSize % 2 == 0) {
		cout << "bad filter-size. pls choose odd-number!" << endl;
		return Mat::zeros(kSize, kSize, CV_32FC1);
	}
	//Mat kernel = Mat::ones(kSize, kSize, CV_32F);
	//for (int x = 0; x < kSize; x++) {
	//	for (int y = 0; y < kSize; y++) {
	//		kernel.col(x).row(y) = kernel.col(x).row(y) / (kSize*kSize);
	//	}
	//}

	//kernel with heigher weight for near pixels//
	int anz = 0;
	int powX = 0;
	int powY = 0;
	Mat kernel = Mat::ones(kSize, kSize, CV_32F);
	for (int x = 0; x < kSize; x++) {
		//proof: before or after anchor
		if (x > ((kSize - 1) / 2)) {
			powX = x - ((kSize - 1) / 2) - 1;
		}
		else {
			powX = x;
		}
		//proof: before or after anchor
		for (int y = 0; y < kSize; y++) {
			if (y > ((kSize - 1) / 2)) {
				powY = y - ((kSize - 1) / 2) - 1;
			}
			else {
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
			kernel.col(x).row(y) = kernel.at<float>(Point(x, y)) / anz;
		}
	}
	return kernel;
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
			//cout << x << ":" << y << " - " << currX << ":" << currY << endl;
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
	////1. circShift Kernel
	////cout << kernel << endl;
	//Mat shiftKernel = circShift(kernel, (kernel.cols - 1) / 2, (kernel.rows - 1) / 2);
	////shiftKernel.row(0).col(0) = 1;
	////shiftKernel.row(0).col(1) = 2;
	////shiftKernel.row(0).col(2) = 3;
	////shiftKernel.row(1).col(0) = 4;
	////shiftKernel.row(1).col(1) = 5;
	////shiftKernel.row(1).col(2) = 6;
	////shiftKernel.row(2).col(0) = 7;
	////shiftKernel.row(2).col(1) = 8;
	////shiftKernel.row(2).col(2) = 9;

	////cout << "shiftKernel: " << endl;
	////cout << shiftKernel << endl ;
	////2. build Kernel in size of image

	//Mat fKernel = Mat::zeros(in.rows, in.cols, CV_32F);
	//repeat(shiftKernel, fKernel.rows / shiftKernel.rows, fKernel.cols / shiftKernel.cols, fKernel);
	//shiftKernel.convertTo(shiftKernel, CV_32F);
	//Mat fImage;// = in.clone();
	//in.convertTo(fImage, CV_32F);
	////cout << "paddedKernel before repeat: " << endl;
	////cout << paddedKernel << endl;
	//
	////cout << "shiftKernel: " << endl;
	////cout << shiftKernel << endl;
	////cout << "fKernel: " << endl;
	////cout << fKernel << endl;
	////cout << "fImage: " << endl;
	////cout << fImage << endl;
	////cout << "paddedKernel after repeat: " << endl;
	////cout << paddedKernel << endl;
	////Mat padded;
	////cout << in << endl;

	////copyMakeBorder(in, padded, 0, m - in.rows, 0, n - in.cols, BORDER_CONSTANT, Scalar::all(0));
	////cout << in.rows << "x" << in.cols << endl;
	////cout << padded << endl;

	////Mat kernelPlanes[] = { Mat_<float>(padded), Mat::zeros(padded.size(), CV_32F) };
	////Mat imagePlanes[] = { Mat_<float>(padded), Mat::zeros(padded.size(), CV_32F) };
	////cout << "kernelPlanes: " << endl;
	////cout << kernelPlanes[0] << endl;
	////cout << "imagePlanes: " << endl;
	////cout << imagePlanes[0] << endl;
	////Mat frequencyKernel;
	////Mat frequencyImage;
	////merge(kernelPlanes, 2, frequencyKernel);         // Add to the expanded another plane with zeros
	////merge(imagePlanes, 2, frequencyImage);
	////cout << "merged frequencyKernel: " << endl;
	////cout << frequencyKernel << endl;

	////dft(frequencyKernel, frequencyKernel);
	////dft(frequencyImage, frequencyImage);
	////// compute the magnitude and switch to logarithmic scale
	////// => log(1 + sqrt(Re(DFT(I))^2 + Im(DFT(I))^2))
	////split(frequencyKernel, kernelPlanes);            // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
	////magnitude(kernelPlanes[0], kernelPlanes[1], kernelPlanes[0]);// planes[0] = magnitude
	////Mat magI = kernelPlanes[0];
	////magI += Scalar::all(1);                    // switch to logarithmic scale
	////log(magI, magI);
	////normalize(magI, magI, 0, 1, CV_MINMAX); // Transform the matrix with float values into a
	////										// viewable image form (float between values 0 and 1).




	////3. dft of kernel and image
	////int m = getOptimalDFTSize(in.rows);
	////int n = getOptimalDFTSize(in.cols); // on the border add zero values

	////Mat imagePlanes[] = { Mat_<float>(fImage), Mat::zeros(fImage.size(), CV_32F) };
	////Mat complexImage;	//Complex plane to contain the DFT coefficients {[0]-Real,[1]-Img}
	////merge(imagePlanes, 2, complexImage);

	////Mat kernelPlanes[] = { Mat_<float>(fKernel), Mat::zeros(fKernel.size(), CV_32F) };
	////Mat complexKernel;	//Complex plane to contain the DFT coefficients {[0]-Real,[1]-Img}
	////merge(kernelPlanes, 2, complexKernel);

	//												  // Reconstructing original imae from the DFT coefficients

	////cout << "paddedKernel: " << endl;
	////cout << paddedKernel << endl;
	//dft(fImage, fImage, 0);
	//dft(fKernel, fKernel,0);
	////cout << "paddedKernel after fourier: " << endl;
	////cout << paddedKernel << endl;

	////4. multiply spectrums of kernel and image in frequenzy domain
	//cout << "fImage: " << endl;
	//cout << fImage << endl;
	//cout << "fKernel: " << endl;
	//cout << fKernel << endl;
	//Mat outputImage;
	//mulSpectrums(fImage, fKernel, outputImage, 0);
	////outputImage = fKernel.mul(fImage);
	//cout << "result mulSpectrums: " << endl;
	//cout << outputImage << endl;
	////outputImage.convertTo(out, CV_8U);
	////5. inverse dft
	////fImage.copyTo(outputImage);
	//dft(outputImage, outputImage, DFT_INVERSE + DFT_SCALE);
	////idft(outputImage, outputImage, DFT_SCALE | DFT_REAL_OUTPUT); // Applying IDFT
	//Mat out;
	//outputImage.convertTo(out, CV_8U);


	////cout << "kernel: " << endl;
	////cout << kernel << endl;
	////cout << "shift kernel: " << endl;
	////cout << shiftKernel << endl;
	////cout << "spectrum kernel: " << endl;
	////cout << fKernel << endl;
	////cout << "spectrum image: " << endl;
	////cout << fImage << endl;
	////cout << "spectrum result: " << endl;
	////cout << outputImage << endl;

	//cout << "spatialConvolution: "<<endl;
	//cout << spatialConvolution(in, shiftKernel)<<endl;
	//cout << "frequencyConvolution: " << endl;
	//cout << out << endl;
return in.clone();//outputImage;
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
			if (fImage.at<float>(Point(x, y)) <= thresh) fImage.row(y).col(x) = 0;
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
Mat Dip3::usm(Mat& in, int type, int size, double thresh, double scale) {

	// some temporary images 
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
		/*Smooth : y0 -> y1  
		2. Subtract : y2 = y0 – y1  
		2.5.Scale : ý2 = s * y2  
		3. Add : y3 = y0 + ý2  
		-> Final(green)*/
	Mat fImage, fIn;
	//original image
	in.convertTo(fIn, CV_32F);
	dft(fIn, fIn, 0);
	//smoothed image
	tmp.convertTo(fImage, CV_32F);
	dft(fImage, fImage, 0);
	//2. subtract

	//cout << "fIn: " << endl;
	//cout << fIn << endl;
	//cout << "fImage: " << endl;
	//cout << fImage << endl;
	subtract(fIn, fImage, fImage);
	//cout << "substract fImage: " << endl;
	//cout << fImage << endl;
	//3. threshold
	fImage = threshold(fImage, thresh);
	//dft(fImage, fImage, DFT_INVERSE + DFT_SCALE);
	//for (int y = 0; y < fImage.rows; y++) {
	//	for (int x = 0; x < fImage.cols; x++) {
	//		if (fImage.at<float>(Point(x, y)) <= thresh) fImage.row(y).col(x) = 0;
	//	}
	//}
	//dft(fImage, fImage, 0);
	
	//4. scale
	fImage = fImage * scale;
	//cout << "scale fImage: " << endl;
	//cout << fImage << endl;
	//5. add
	//fImage += in;
	add(fImage, fIn, fImage);
	//cout << "add fImage: " << endl;
	//cout << fImage << endl;
	//cin.get();
	Mat outputImage = in.clone();// Mat::zeros(fImage.rows, fImage.cols, CV_8U);
	dft(fImage, fImage, DFT_INVERSE + DFT_SCALE);
	//idft(outputImage, outputImage, DFT_SCALE | DFT_REAL_OUTPUT); // Applying IDFT

	//fImage.convertTo(outputImage, CV_8U);
	return fImage;

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
	Mat mat1 = Mat(3, 3, CV_32F);
	Mat mat2 = Mat(3, 3, CV_32F);
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
			}
			else {
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
	float ref[9][9] = { {0, 0, 0, 0, 0, 0, 0, 0, 0},
					   {0, 1, 1, 1, 1, 1, 1, 1, 0},
					   {0, 1, 1, 1, 1, 1, 1, 1, 0},
					   {0, 1, 1, (8 + 255) / 9., (8 + 255) / 9., (8 + 255) / 9., 1, 1, 0},
					   {0, 1, 1, (8 + 255) / 9., (8 + 255) / 9., (8 + 255) / 9., 1, 1, 0},
					   {0, 1, 1, (8 + 255) / 9., (8 + 255) / 9., (8 + 255) / 9., 1, 1, 0},
					   {0, 1, 1, 1, 1, 1, 1, 1, 0},
					   {0, 1, 1, 1, 1, 1, 1, 1, 0},
					   {0, 0, 0, 0, 0, 0, 0, 0, 0} };
	for (int y = 1; y < 8; y++) {
		for (int x = 1; x < 8; x++) {
			if (abs(output.at<float>(y, x) - ref[y][x]) > 0.0001) {
				cout << "ERROR: Dip3::frequencyConvolution(): Convolution result contains wrong values!" << endl;
				return;
			}
		}
	}
	cout << "Message: Dip3::frequencyConvolution() seems to be correct" << endl;
}
