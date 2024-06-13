#ifndef MG_IMAGE_H_
#define MG_IMAGE_H_

#include<iostream>
#include<opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <vector>

bool read_image(const char* filename, cv::Mat& img, int flag=1);
void get_channel_blue(cv::Mat inputImg, cv::Mat& outputImg);
void get_channel_green(cv::Mat inputImg, cv::Mat& outputImg);
void get_channel_red(cv::Mat inputImg, cv::Mat& outputImg);
void get_channels(cv::Mat inputImg, cv::Mat* channelImg);
void get_channel_brightness(cv::Mat input_hsv_img, cv::Mat& outImg);
void image_pyramid_down(cv::Mat inputImg, std::vector<cv::Mat>& pyramid, int depth, float scale);
void auto_adjust_rgb_color(cv::Mat inputImg, cv::Mat& outputImg, float k, float shift);
void auto_adjust_rgb_color(cv::Mat inputImg, cv::Mat& outputImg);
void auto_adjust_gray(cv::Mat inputImg, cv::Mat& outputImg);

void warp_triangle(cv::Mat inputImg, cv::Mat& outputImg, std::vector<cv::Point2f> inTri, std::vector<cv::Point2f> outTri);
void warp_triangle_mask(cv::Mat inputImg, cv::Rect& bbox_out, cv::Mat& mask, cv::Mat& warped_triangle, std::vector<cv::Point2f> inTri, std::vector<cv::Point2f> outTri);

#endif