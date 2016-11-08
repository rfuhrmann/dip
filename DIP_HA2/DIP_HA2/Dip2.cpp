//============================================================================
// Name        : Dip2.cpp
// Author      : Ronny Haensch
// Version     : 2.0
// Copyright   : -
// Description : 
//============================================================================

#include "Dip2.h"
#include "list"

// convolution in spatial domain
/*
src:     input image
kernel:  filter kernel
return:  convolution result
*/
Mat Dip2::spatialConvolution(Mat& src, Mat& kernel){
   // TO DO !!
	//cout << src << endl;
	//cout << endl << kernel << endl;
	//1. flip kernel
	Mat flipKernel = Mat(kernel.cols, kernel.rows, CV_32FC1);
	for (int x = 0; x < (kernel.cols); x++) {
		//for (int x = 0; x < (mat.cols - 1) / 2; x++) {
		for (int y = 0; y < (kernel.rows); y++) {
			//for (int y = 0; y < (mat.rows - 1) / 2; y++) {
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
				//cout << mean(mat2).val[0] << endl;
				//outputImage.at<float>(Point(x, y)) = mean(mat2).val[0]*(kernel.cols*kernel.rows);
				outputImage.col(x).row(y) = mean(mat2)*(kernel.cols*kernel.rows);
			};
		}
	}

	////Convolution-Method of OpenCV
	//int ddepth = -1;
	//Point anchor = Point(-1, -1);
	//double delta = 0;
	//filter2D(src, outputImage, ddepth, flipKernel, anchor, delta, BORDER_DEFAULT);
	//imshow("filter2D", outputImage);
	//cout << outputImage << endl;
   return outputImage;
}

Point2i getMedian(Mat& src, Point2i anchor, int size) {
	//cout <<"Anchor: "<< anchor.x << "x" << anchor.y << endl;
	// proof for valid input
	//if (src.empty() || anchor.x == NULL || anchor.y == NULL || size == NULL) return -1;
	//if (size % 2 == 0) {
	//	cout << "Pls enter a odd-numbered size to enable anchor!" << endl;
	//	return -1;
	//}
	Rect rect(anchor.x - ((size - 1)) / 2, anchor.y-((size-1)/2), size, size);
	Mat mat = src(rect);


	Point2i medianPixel = anchor;
	Mat kernel = Mat::zeros(size, size, CV_8UC2);
	int srcX = 0;
	int srcY = 0;
	//float array[9];
	list<float> list;
	/*cout << "List1: " << list.size << endl;
	list.push_back(2);
	cout << "List2: " << list.size << endl;
	*///for (int x = anchor.x - ((size - 1) / 2); x <= anchor.x + ((size - 1) / 2); x++) {
	for (int x = 0; x < size; x++) {
		//for (int y = anchor.y - ((size - 1) / 2); y <= anchor.y + ((size - 1) / 2); y++) {
		if (anchor.x - ((size - 1) / 2) + x < 0) {
			srcX = 0;
		}
		else if (anchor.x - ((size - 1) / 2) + x >= src.cols) {
			srcX = src.cols - 1;
		}
		else {
			srcX = anchor.x - ((size - 1) / 2) + x;
		}
		for (int y = 0; y < size; y++) {

			if (anchor.y - ((size - 1) / 2) + y < 0) {
				srcY = 0;
			}
			else if (anchor.y - ((size - 1) / 2) + y >= src.rows) {
				srcY = src.rows - 1;
			}
			else {
				srcY = anchor.y - ((size - 1) / 2) + y;
			}
			list.push_back(src.at<float>(Point(srcX, srcY)));
			//.push_back(Point(srcX, srcY));
		}
	}
	list.sort();
	//cout << "list: " << integerList.size << endl;
	for (int i = 0; i < size; i++) {
		list.pop_back();
	}
	/*while (src.at<float>(point2iList.back()) != list.back()) {
		if (point2iList.empty()) {
			cout << "ERROR - point2iList is empty" << endl;
			waitKey(0);
			return anchor;
		}
		point2iList.pop_back();
	}*/
	//cout << "medien: " << point2iList.back() << endl;

	return Point(1, 1);//point2iList.back();
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
	//cout << "average:" << endl;
	//cout << src << endl;
	Mat kernel = Mat::ones(kSize, kSize, CV_32F);
	for (int x = 0; x < kSize; x++) {
		for (int y = 0; y < kSize; y++) {
			kernel.col(x).row(y) = kernel.col(x).row(y) / (kSize*kSize);
		}
	}
	Mat outputImage = spatialConvolution(src, kernel);
	//cout << outputImage << endl;
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
	//Mat ownImage = imread("noiseType_1.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	Mat outputImage = src.clone();//ownImage.clone();
	Point2i anchor;
	//cout << "kSize: " << kSize << endl;
	//cout << "srcSize: " << src.cols << "x" << src.rows << endl;
	for (int x = 0; x < src.cols; x++) {
		anchor.x = x;
		for (int y = 0; y < src.rows; y++) {
			anchor.y = y;
			//outputImage.at<float>(Point(x, y)) = src.at<float>(getMedian(src, anchor, kSize));
		}
	}
	//namedWindow("Test");
	//imshow("Test", outputImage);
	Mat dst;
	medianBlur(src, dst, kSize);
	imshow("source", src);
	imshow("result", dst);
	return dst;//src.clone();
}

// the bilateral filter
/*
src:     input image
kSize:   size of the kernel --> used to compute std-dev of spatial kernel
sigma:   standard-deviation of the radiometric kernel
return:  filtered image
*/
Mat Dip2::bilateralFilter(Mat& src, int kSize, double sigma){
  
    return src.clone();

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
	Mat restorated1 = noiseReduction(noise1, "median", 3);
	Mat restorated2 = noiseReduction(noise2, "", 1);
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
