#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

int invalid = -999999;
// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    float dobsigmasq = 2 * pow(sigma , 2);
    float denom = sqrt(M_PI * dobsigmasq);
    int w = 2 * ( (int) ceil(3 * sigma)) + 1;
    image filter = make_image(1,w,1);
    int offset = w / 2;
    for(int i = 0; i < w; i++){
        int x = i - offset;
        float exparg = - (pow(x, 2)) / dobsigmasq;
        float val = exp(exparg) / denom;
        set_pixel(filter, 0, i, 0, val);
    }
    return filter;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    if(0){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        // TODO: optional, use two convolutions with 1d gaussian filter.
        // If you implement, disable the above if check.
        image f1d = make_1d_gaussian(sigma);
        image horizsmoth = convolve_image(im, f1d, 1);
        f1d.w = f1d.h;
        f1d.h= 1;
        image s = convolve_image(horizsmoth, f1d, 1);
        free_image(f1d);
        free_image(horizsmoth);
        return s;
    }
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image gs = make_image(im.w, im.h, 3);
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image gx = convolve_image(im, gx_filter, 0);
    /*feature_normalize(gx);*/
    image gy = convolve_image(im, gy_filter, 0);
    /*feature_normalize(gy);*/
    float x;
    float y;
    for (int i = 0; i < im.w; i++){
      for(int j = 0; j < im.h; j++){
        x = get_pixel(gx, i, j, 0);
        y = get_pixel(gy, i, j, 0);
        set_pixel(gs, i, j, 0, x * x);
        set_pixel(gs, i, j, 1, y * y);
        set_pixel(gs, i, j, 2, x * y);
      }
    }
    /*feature_normalize(gs);*/
    /*image gaussian = make_gaussian_filter(sigma);*/
    /*image S = convolve_image(gs, gaussian, 1);*/
    image S = smooth_image(gs, sigma);
    free_image(gs);
    free_image(gx_filter);
    free_image(gy_filter);
    free_image(gx);
    free_image(gy);

    // TODO: calculate structure matrix for im.
    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    float xx;
    float xy;
    float yy;
    float val;
    float det;
    float trace;
    for (int i = 0; i < S.w; i++){
      for (int j = 0; j < S.h; j++){
        xx = get_pixel(S, i, j, 0);
        yy = get_pixel(S, i, j, 1);
        xy = get_pixel(S, i, j, 2);
        trace = xx + yy;
        det = xx * yy - xy * xy;
        val = det - 0.06 * trace * trace;
        set_pixel(R, i, j, 0, val);
      }
    }

    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
  image r = copy_image(im);
  // TODO: perform NMS on the response map.
  // for every pixel in the image:
  //     for neighbors within w:
  //         if neighbor response greater than pixel response:
  //             set response to be very low (I use -999999 [why not 0??])
  int set_invalid;
  int x, y;
  float curr;
  float walker;
  float elmas = -99999999.0;
  for (int i = 0; i < im.w; i++){
    for (int j = 0; j < im.h; j++){
      set_invalid = 0;
      curr = get_pixel(im, i, j, 0);
      if (curr> elmas) {
        elmas = curr;
      }
      for (int k = -w; (k <=w) & (set_invalid == 0); k++){
        for (int l = -w; (l <=w) & (set_invalid == 0); l++){
          x = i + k;
          y = j + l;
          if ((x>=0) & (x < im.w) & (y >= 0) & (y < im.h)){
            walker = get_pixel(im, x, y, 0);
            if (walker > curr) {
              set_pixel(r, i, j, 0, invalid);
              set_invalid = 1;
            }
          }
        }
      }
      if (set_invalid == 0) {
        set_invalid = 0;
      }
    }
  }
  return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);


    //TODO: count number of responses over threshold
    int count = 0; // change this
    for (int i = 0; i < (Rnms.w * Rnms.h * Rnms.c); i++)
    {
      if (Rnms.data[i] > thresh) {
        count++;
      }
    }


    
    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.
    int j = 0;
    for (int i = 0; i < (Rnms.w * Rnms.h * Rnms.c); i++){
      if (Rnms.data[i] > thresh) {
      d[j] = describe_index(im, i);
      j++;
      }}


    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
