#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // TODO
    float total = 0;
    for(int i=0; i < im.w; i++){
      for(int j = 0; j < im.h; j++){
        for(int c = 0; c < im.c; c++){
        /*total += im.data[index_calc(im.w, im.h, i, j, c)];*/
        total += get_pixel(im, i, j, c);
        }
      }
    }
    float val;
    for(int i=0; i < im.w; i++){
      for(int j = 0; j < im.h; j++){
        for(int c = 0; c < im.c; c++){
        /*val = im.data[index_calc(im.w, im.h, i, j, c)];*/
        val = get_pixel(im, i, j, c) / total;
        /*im.data[index_calc(im.w, im.h, i, j, c)] = val / total;*/
        set_pixel(im, i, j, c, val);
        }
      }
    }
}

image make_box_filter(int w)
{
    // TODO
    image im = make_image(w,w,1);
    for (int i = 0; i < w; i++) {
      for (int j = 0; j < w; j++) {
        set_pixel(im, i, j, 0, 1.0);
      }
    }
    l1_normalize(im);
    return im;
}

image convolve_image(image im, image filter, int preserve)
{
    assert((filter.c == 1) | (filter.c == im.c));
    image new_im;
    int i_offset;
    int j_offset;
if ( (im.c == filter.c) & (preserve == 0) ) {
      new_im = make_image(im.w, im.h, 1);
      for(int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
          float val = 0;
          for (int c = 0; c < im.c; c++) {
            i_offset = filter.w / 2;
            j_offset = filter.h / 2;
            float imgval;
            float filtval;
            for (int k = 0; k < filter.w; k++) {
              for (int l = 0; l < filter.h; l++){
                imgval = get_pixel(im, i - i_offset + k, j - j_offset + l, c);
                filtval = get_pixel(filter,k, l, c);
                val += imgval * filtval;
              }
            }
          }
          set_pixel(new_im, i, j, 0, val);
        }
      }
    } else if ((filter.c == 1) & (preserve == 0)) {
      new_im = make_image(im.w, im.h, 1);
      for(int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float val = 0;
          for (int c = 0; c < im.c; c++) {
            i_offset = filter.w / 2;
            j_offset = filter.h / 2;
            float imgval;
            float filtval;
            for (int k = 0; k < filter.w; k++) {
              for (int l = 0; l < filter.h; l++){
                imgval = get_pixel(im, i - i_offset + k, j - j_offset + l, c);
                filtval = get_pixel(filter,k, l, 0);
                val += imgval * filtval;
              }
            }
          }
            /*val = fmax(val, 0.0);*/
            /*val = fmin(val, 1.0);*/
            set_pixel(new_im, i, j, 0, val);
            /*set_pixel(new_im, i, j, 1, val);*/
            /*set_pixel(new_im, i, j, 2, val);*/
        }
      }
    } else if ((filter.c == 1) & (preserve == 1)) {
      new_im = make_image(im.w, im.h, im.c);
      for(int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
          for (int c = 0; c < im.c; c++) {
            i_offset = filter.w / 2;
            j_offset = filter.h / 2;
            float val = 0;
            float imgval;
            float filtval;
            for (int k = 0; k < filter.w; k++) {
              for (int l = 0; l < filter.h; l++){
                imgval = get_pixel(im, i - i_offset + k, j - j_offset + l, c);
                filtval = get_pixel(filter,k, l, 0);
                val += imgval * filtval;
              }
            }
            /*val = fmax(val, 0.0);*/
            /*val = fmin(val, 1.0);*/
            set_pixel(new_im, i, j, c, val);
          }
        }
      }
    } else {
      new_im = make_image(im.w, im.h, im.c);
      printf("ERROOOOR");
    }
    return new_im;


    // TODO
    /*return make_image(1,1,1);*/
}

image make_highpass_filter()
{
    // TODO
    image im = make_image(3,3,1);
    set_pixel(im, 0, 0, 0, 0.0);
    set_pixel(im, 0, 1, 0, -1.0);
    set_pixel(im, 0, 2, 0, 0.0);
    set_pixel(im, 1, 0, 0, -1.0);
    set_pixel(im, 1, 1, 0, 4.0);
    set_pixel(im, 1, 2, 0, -1.0);
    set_pixel(im, 2, 0, 0, 0.0);
    set_pixel(im, 2, 1, 0, -1.0);
    set_pixel(im, 2, 2, 0, 0.0);
    /*set_pixel(im, 0, 0, 0, 0.0);*/
    /*set_pixel(im, 0, 1, 0, 0.0);*/
    /*set_pixel(im, 0, 2, 0, 0.0);*/
    /*set_pixel(im, 1, 0, 0, 0.0);*/
    /*set_pixel(im, 1, 1, 0, 0.0);*/
    /*set_pixel(im, 1, 2, 0, 0.0);*/
    /*set_pixel(im, 2, 0, 0, 0.0);*/
    /*set_pixel(im, 2, 1, 0, 1.0);*/
    /*set_pixel(im, 2, 2, 0, 0.0);*/
    return im;
}

image make_sharpen_filter()
{
    image im = make_image(3,3,1);
    set_pixel(im, 0, 0, 0, 0.0);
    set_pixel(im, 0, 1, 0, -1.0);
    set_pixel(im, 0, 2, 0, 0.0);
    set_pixel(im, 1, 0, 0, -1.0);
    set_pixel(im, 1, 1, 0, 5.0);
    set_pixel(im, 1, 2, 0, -1.0);
    set_pixel(im, 2, 0, 0, 0.0);
    set_pixel(im, 2, 1, 0, -1.0);
    set_pixel(im, 2, 2, 0, 0.0);
    return im;
}

image make_emboss_filter()
{
    image im = make_image(3,3,1);
    set_pixel(im, 0, 0, 0, -2.0);
    set_pixel(im, 0, 1, 0, -1.0);
    set_pixel(im, 0, 2, 0, 0.0);
    set_pixel(im, 1, 0, 0, -1.0);
    set_pixel(im, 1, 1, 0, 1.0);
    set_pixel(im, 1, 2, 0, 1.0);
    set_pixel(im, 2, 0, 0, 0.0);
    set_pixel(im, 2, 1, 0, 1.0);
    set_pixel(im, 2, 2, 0, 2.0);
    return im;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: Preserve is for emboss and sharpen, highpass shouldn't preserve the channels

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: cap values lesser than 0 or greater than 1.

image make_gaussian_filter(float sigma)
{
    // TODO
    float dobsigmasq = 2 * pow(sigma , 2);
    float denom = M_PI * dobsigmasq;
    int w = 2 * ( (int) ceil(3 * sigma)) + 1;
    image filter = make_image(w,w,1);
    int offset = w / 2;
    for(int i = 0; i < w; i++){
      for (int j = 0; j < w; j++){
        int x = i - offset;
        int y = j - offset;
        float exparg = - (pow(x, 2) + pow(y, 2)) / dobsigmasq;
        float val = exp(exparg) / denom;
        set_pixel(filter, i, j, 0, val);
      }
    }
    return filter;
}

image add_image(image a, image b)
{
    assert(a.w == b.w);
    assert(a.h == b.h);
    assert(a.c == b.c);
    image im = make_image(a.w, a.h, a.c);
    float aval;
    float bval;
    float val;
    for (int i = 0; i < a.w; i++) {
      for (int j = 0; j < a.h; j++) {
        for (int c = 0; c < a.c; c++) {
          aval = get_pixel(a, i, j, c);
          bval = get_pixel(b, i, j, c);
          val = aval + bval;
          val = fmin(val, 1.0);
          val = fmax(val, 0.0);
          set_pixel(im, i, j, c, val);
        }
      }
    }
    return im;
}

image sub_image(image a, image b)
{
    assert(a.w == b.w);
    assert(a.h == b.h);
    assert(a.c == b.c);
    image im = make_image(a.w, a.h, a.c);
    float aval;
    float bval;
    float val;
    for (int i = 0; i < a.w; i++) {
      for (int j = 0; j < a.h; j++) {
        for (int c = 0; c < a.c; c++) {
          aval = get_pixel(a, i, j, c);
          bval = get_pixel(b, i, j, c);
          val = aval - bval;
          /*val = fmin(val, 1.0);*/
          /*val = fmax(val, 0.0);*/
          set_pixel(im, i, j, c, val);
        }
      }
    }
    return im;
}

image make_gy_filter()
{
    image filter = make_image(3,3,1);
    set_pixel(filter, 0, 0, 0, -1.0);
    set_pixel(filter, 0, 1, 0, 0.0);
    set_pixel(filter, 0, 2, 0, 1.0);
    set_pixel(filter, 1, 0, 0, -2.0);
    set_pixel(filter, 1, 1, 0, 0.0);
    set_pixel(filter, 1, 2, 0, 2.0);
    set_pixel(filter, 2, 0, 0, -1.0);
    set_pixel(filter, 2, 1, 0, 0.0);
    set_pixel(filter, 2, 2, 0, 1.0);
    return filter;
}

image make_gx_filter()
{
    image filter = make_image(3,3,1);
    set_pixel(filter, 0, 0, 0, -1.0);
    set_pixel(filter, 0, 1, 0, -2.0);
    set_pixel(filter, 0, 2, 0, -1.0);
    set_pixel(filter, 1, 0, 0, 0.0);
    set_pixel(filter, 1, 1, 0, 0.0);
    set_pixel(filter, 1, 2, 0, 0.0);
    set_pixel(filter, 2, 0, 0, 1.0);
    set_pixel(filter, 2, 1, 0, 2.0);
    set_pixel(filter, 2, 2, 0, 1.0);
    return filter;
}

void feature_normalize(image im)
{
  float minval = 1000;
  float maxval = -1111;
  float pixval;
  for(int i = 0; i < im.w; i++) {
    for (int j = 0; j < im.h; j++) {
      for (int c = 0; c < im.c; c++) {
        pixval = get_pixel(im, i, j, c);
        minval = fmin(minval, pixval);
        maxval = fmax(maxval, pixval);
      }
    }
  }
  float range = maxval - minval;
  for(int i = 0; i < im.w; i++) {
    for (int j = 0; j < im.h; j++) {
      for (int c = 0; c < im.c; c++) {
        if (range == 0) {
          set_pixel(im, i, j, c, 0);
        } else {
          pixval = get_pixel(im, i, j, c);
          pixval -= minval;
          pixval /= range;
          pixval /= im.c;
          set_pixel(im, i, j, c, pixval);
        }
      }
    }
  }
}

image *sobel_image(image im)
{
    image *impo = calloc(2, sizeof(image));
    impo[0].data = calloc(im.w * im.h, sizeof(float));
    impo[1].data = calloc(im.w * im.h, sizeof(float));
    impo[0].w = im.w;
    impo[0].h = im.h;
    impo[0].c = 1;
    impo[1].w = im.w;
    impo[1].h = im.h;
    impo[1].c = 1;
    image gx = make_gx_filter();
    image gy = make_gy_filter();
    image x = convolve_image(im, gx, 0);
    image y = convolve_image(im, gy, 0);
    float xval, yval, magn, angle;
    for (int i = 0; i < im.w; i++) {
      for (int j = 0; j < im.h; j++) {
        xval = get_pixel(x, i, j, 0);
        yval = get_pixel(y, i, j, 0);
        magn = sqrt(pow(xval, 2) + pow(yval, 2));
        /*magn = fmin(fmax(magn, 0.0), 1.0);*/
        set_pixel(impo[0], i, j, 0, magn);
        angle = atan2(yval, xval);
        /*angle = angle * 6 / ( 2 * M_PI);*/
        /*angle = fmin(fmax(magn, 0.0), 6.0);*/
        set_pixel(impo[1], i, j, 0, angle);
      }
    }
    free_image(gx);
    free_image(gy);
    free_image(x);
    free_image(y);
    /*feature_normalize(x);*/
    /*feature_normalize(y);*/
    /*save_image(x, "ladegx");*/
    /*save_image(y, "ladegy");*/
    return impo;
}

image colorize_sobel(image im)
{
    image *sobel = sobel_image(im);
    image res = make_image(im.w, im.h, 3);
    feature_normalize(sobel[0]);
    feature_normalize(sobel[1]);
    float magn;
    float angle;
    for (int i = 0; i < im.w; i++) {
      for (int j = 0; j < im.h; j++) {
        magn = get_pixel(sobel[0], i, j, 0);
        angle = get_pixel(sobel[1], i, j, 0);
        set_pixel(res, i, j, 0, angle);
        set_pixel(res, i, j, 1, magn);
        set_pixel(res, i, j, 2, magn);
        }
      }
    hsv_to_rgb(res);
    free_image(*sobel);
    return res;
}
