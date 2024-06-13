#include "image.h"

#include <opencv2/opencv.hpp>

void warp_triangle(cv::Mat inputImg, cv::Mat& outputImg, std::vector<cv::Point2f> inTri, std::vector<cv::Point2f> outTri)
{
    // Check if the input image is a color image with floating numbers(0<= intensity <=1)
    if (inputImg.type() != CV_32FC3)
    {
        std::cout << "Error, the type of input image is not supported for this triangle warping method" << std::endl;
        return;
    }

    // Set the output image as the white first
    //cv::Mat tempOutputImg = cv::Mat(inputImg.size(), inputImg.type());
    //tempOutputImg = cv::Scalar(1.0, 1.0, 1.0);

    // Find the bounding boxes for triangles
    cv::Rect bbox1 = cv::boundingRect(inTri);
    cv::Rect bbox2 = cv::boundingRect(outTri);

    // Triangles in the local coordinates of the bounding boxes
    std::vector<cv::Point2f> inTriLocal(3), outTriLocal(3);
    for (int i = 0; i < 3; ++i)
    {
        inTriLocal[i] = cv::Point2f(inTri[i].x - bbox1.x, inTri[i].y - bbox1.y);
        outTriLocal[i] = cv::Point2f(outTri[i].x - bbox2.x, outTri[i].y - bbox2.y);
    }

    // Cropped image by the bounding box in the input image
    cv::Mat inBBoxImg;
    inputImg(bbox1).copyTo(inBBoxImg);

    // Affine transform from input triangle to output triangle
    cv::Mat A = cv::getAffineTransform(inTriLocal, outTriLocal);

    // Apply the affine transform to obtain the bounding box image in the output image
    cv::Mat outBBoxImg = cv::Mat::zeros(bbox2.height, bbox2.width, inBBoxImg.type());
    cv::warpAffine(inBBoxImg, outBBoxImg, A, outBBoxImg.size(), cv::INTER_LINEAR, cv::BORDER_REFLECT_101);

    // Find the mask for the triangle
    cv::Mat mask = cv::Mat::zeros(outBBoxImg.size(), outBBoxImg.type());
    std::vector<cv::Point> outTriLocalInt(3);
    for (int i = 0; i < 3; ++i)
    {
        outTriLocalInt[i] = cv::Point((int)outTriLocal[i].x, (int)outTriLocal[i].y);
    }

    cv::fillConvexPoly(mask, outTriLocalInt, cv::Scalar(1.0, 1.0, 1.0), 16, 0);

    // Masked triangle image
    cv::multiply(outBBoxImg, mask, outBBoxImg);
    cv::multiply(outputImg(bbox2), cv::Scalar(1.0, 1.0, 1.0) - mask, outputImg(bbox2));
    outputImg(bbox2) = outputImg(bbox2) + outBBoxImg;

    // For visualization
    /*cv::Mat tempImg;
    outputImg.convertTo(tempImg, CV_8UC3, 255.9);
    cv::imwrite("warp_triangle.jpg", tempImg);*/
}

void warp_triangle_mask(cv::Mat inputImg, cv::Rect& bbox_out, cv::Mat& mask, cv::Mat& warped_triangle, std::vector<cv::Point2f> inTri, std::vector<cv::Point2f> outTri)
{
    // Check if the input image is a color image with floating numbers(0<= intensity <=1)
    if (inputImg.type() != CV_32FC3 )
    {
        std::cout << "Error, the type of input image is not supported for this triangle warping method" << std::endl;
        return;
    }

    // Set the output image as the white first
    //cv::Mat tempOutputImg = cv::Mat(inputImg.size(), inputImg.type());
    //tempOutputImg = cv::Scalar(1.0, 1.0, 1.0);

    // Find the bounding boxes for triangles
    cv::Rect bbox1 = cv::boundingRect(inTri);
    cv::Rect bbox2 = cv::boundingRect(outTri);
    bbox_out = bbox2;

    // Triangles in the local coordinates of the bounding boxes
    std::vector<cv::Point2f> inTriLocal(3), outTriLocal(3);
    for (int i = 0; i < 3; ++i)
    {
        inTriLocal[i] = cv::Point2f(inTri[i].x - bbox1.x, inTri[i].y - bbox1.y);
        outTriLocal[i] = cv::Point2f(outTri[i].x - bbox2.x, outTri[i].y - bbox2.y);
    }

    // Cropped image by the bounding box in the input image
    cv::Mat inBBoxImg;
    inputImg(bbox1).copyTo(inBBoxImg);

    // Affine transform from input triangle to output triangle
    cv::Mat A = cv::getAffineTransform(inTriLocal, outTriLocal);

    // Apply the affine transform to obtain the bounding box image in the output image
    cv::Mat outBBoxImg = cv::Mat::zeros(bbox2.height, bbox2.width, inBBoxImg.type());
    cv::warpAffine(inBBoxImg, outBBoxImg, A, outBBoxImg.size(), cv::INTER_LINEAR, cv::BORDER_REFLECT_101);

    // Find the mask for the triangle
    mask = cv::Mat::zeros(outBBoxImg.size(), outBBoxImg.type());
    std::vector<cv::Point> outTriLocalInt(3);
    for (int i = 0; i < 3; ++i)
    {
        outTriLocalInt[i] = cv::Point((int)outTriLocal[i].x, (int)outTriLocal[i].y);
    }

    cv::fillConvexPoly(mask, outTriLocalInt, cv::Scalar(1.0, 1.0, 1.0), 16, 0);

    // Masked triangle image
    cv::multiply(outBBoxImg, mask, outBBoxImg);
    outBBoxImg.copyTo(warped_triangle);
    //cv::multiply(warped_triangle, cv::Scalar(1.0, 1.0, 1.0) - mask, outputImg(bbox2));
    //outputImg(bbox2) = outputImg(bbox2) + outBBoxImg;

    //// For visualization
    //cv::Mat tempImg;
    //outBBoxImg.convertTo(tempImg, CV_8UC3, 255.9);
    //cv::imwrite("warp_triangle.jpg", tempImg);
}

bool read_image(const char* filename, cv::Mat& img, int flag)
{
    std::cout << "- Reading image : " << filename;

    img = cv::imread(filename, flag);

    if (img.empty())
    {
        std::cout << "...ERROR" << std::endl;
        return false;
    }

    std::cout << std::endl;

    return true;
}

void image_pyramid_down(cv::Mat inputImg, std::vector<cv::Mat>& pyramid, int depth, float scale)
{
    pyramid.resize(depth + 1);

    inputImg.copyTo(pyramid[0]);

    for (int i = 0; i < depth; ++i)
    {
        cv::pyrDown(pyramid[i], pyramid[i + 1], cv::Size((int)(pyramid[i].cols*scale), (int)(pyramid[i].rows*scale)));
    }
}

void get_channel_blue(cv::Mat inputImg, cv::Mat& outputImg)
{
    cv::Mat ch[3];
    get_channels(inputImg, ch);
    ch[0].copyTo(outputImg);
}
void get_channel_green(cv::Mat inputImg, cv::Mat& outputImg)
{
    cv::Mat ch[3];
    get_channels(inputImg, ch);
    ch[1].copyTo(outputImg);
}
void get_channel_red(cv::Mat inputImg, cv::Mat& outputImg)
{
    cv::Mat ch[3];
    get_channels(inputImg, ch);
    ch[2].copyTo(outputImg);
}
void get_channels(cv::Mat inputImg, cv::Mat* channelImg)
{
    cv::split(inputImg, channelImg);
}

void get_channel_brightness(cv::Mat inputImg, cv::Mat& out_brightness_channel)
{
    cv::Mat temp;
	cv::cvtColor(inputImg, temp, cv::COLOR_RGB2HSV);
    cv::Mat ch[3];
    get_channels(temp, ch);
    ch[2].copyTo(out_brightness_channel);
}

void auto_adjust_rgb_color(cv::Mat inputImg, cv::Mat& outputImg, float k, float shift)
{
    const int w = inputImg.cols;
    const int h = inputImg.rows;

    inputImg.copyTo(outputImg);

    //int max_intensity = 0;

    int ch = outputImg.channels();

    std::vector<float> max_intensity(ch, 0.0f), mean(ch), dev(ch);

    unsigned char* data = outputImg.data;
    for (int i = 0; i < w*h; ++i)
    {
        for (int j = 0; j < ch; ++j)
        {
            unsigned char pix_value = data[ch*i + j];
            if (pix_value > max_intensity[j])   max_intensity[j] = pix_value;

            mean[j] += pix_value;
            dev[j] += pix_value*pix_value;
        }
    }

    for (int j = 0; j < ch; ++j)
    {
        mean[j] /= w*h;
        dev[j] = sqrt(dev[j] / (w*h) - mean[j] * mean[j]);
        std::cout << "ch[" << j << "]: mean=" << mean[j] << ", dev=" << dev[j] << std::endl;
    }

    std::vector<float> cutoff(ch);
    float max_cutoff = 0.0f;
    int max_ch = 0;
    for (int j = 0; j < ch; ++j)
    {
        cutoff[j] = mean[j] + k*dev[j];
        if (cutoff[j] > 255.0)  cutoff[j] = 255.0;
        if (cutoff[j] > max_cutoff)
        {
            max_cutoff = cutoff[j];
            max_ch = j;
        }
    }

    std::cout << "max_cutoff=" << max_cutoff << std::endl;

    /*std::vector<cv::Mat> chImg(ch), outChImg(ch);
    cv::split(inputImg, chImg);

    for (int j = 0; j < ch; ++j)
    {
        chImg[j].convertTo(outChImg[j], -1, 255 / cutoff[j], 0.0);
    }

    cv::merge(outChImg, outputImg);//*/

    inputImg.convertTo(outputImg, -1, 255 / max_cutoff, shift);

    /*std::vector<int> count(ch, 0);
    for (int i = 0; i < w*h; ++i)
    {
        for (int j = 0; j < ch; ++j)
        {
            float pix_value = (float) data[ch*i + j];
            if (pix_value >= max_cutoff)
            {
                data[ch*i + j] = 255;
                count[j] ++;
            }
            else
            //if (pix_value <= max_cutoff)
            {
                //data[ch*i + j] = (int)(pix_value / max_cutoff * 225.9999f);
                data[ch*i + j] = (int)(pix_value / max_cutoff * 225.9999f);
            }
            //else if (j == max_ch)
            //{
            //    count++;
            //}
        }// end for(j)
    }// end for(i)

    const int numpixels = w*h;
    std::cout << "max_ch =" << max_ch << ", count=" << count[0] << "(" << count[0] / (float)numpixels *100<< ")"
        << ", " << count[1] << "(" << count[1] / (float)numpixels * 100 << ")" << ", " 
        << count[2] << "(" << count[2] / (float)numpixels * 100 << ")" << std::endl;
        */
}

void auto_adjust_gray(cv::Mat inputImg, cv::Mat& outputImg)
{
    std::cout << "adjusting...." << std::endl;

    const int w = inputImg.cols;
    const int h = inputImg.rows;

    const int num_pixels = w*h;

    // Compute histogram
    std::vector<int> histo(256, 0);
    unsigned char* data = inputImg.data;

    for (int i = 0; i < num_pixels; ++i)
    {
        //std::cout << (int)data[i] << " ";
        histo[data[i]]++;
    }

    // Find the cuts (1%)
    // upper cut
    std::cout << "Finding cuts..." << std::endl;
    int count[2], cut[2];
    count[0] = 0;
    for (int i = 0; i < 255; ++i)
    {
        count[0] += histo[255 - i];
        if (count[0] * 100.0f / num_pixels > 0.5f)
        {
            cut[0] = 255 - i;
            break;
        }
    }

    // lower cut
    count[1] = 0;
    for (int i = 0; i < 255; ++i)
    {
        count[1] += histo[i];
        if (count[1] * 100.0f / num_pixels > 0.5f)
        {
            cut[1] = i;
            break;
        }
    }

    // Find the transform factors
    float alpha = 255.0f / (cut[0] - cut[1]);
    float beta = -255.0f*cut[1] / (cut[0] - cut[1]);

    std::cout << "cut : " << cut[1] << "~" << cut[0] << ", alpha=" << alpha << ", beta=" << beta << std::endl;

    inputImg.convertTo(outputImg, -1, alpha, beta);
}

void auto_adjust_rgb_color(cv::Mat inputImg, cv::Mat& outputImg)
{
    int num_channels = inputImg.channels();
    std::vector<cv::Mat> ch_img(num_channels), out_ch_img(num_channels);
    cv::split(inputImg, ch_img);

    for (int i = 0; i < num_channels; ++i)
    {
        auto_adjust_gray(ch_img[i], out_ch_img[i]);
    }

    cv::merge(out_ch_img, outputImg);
}