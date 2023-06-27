#include <math.h>
#include <stdio.h>
#include "image.h"

int index_calc(int w,int  h,int  i,int  j,int  c){
  return c * (w * h) + w * j + i;
}

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO Fill in
    int xi = round(x);
    int yi = round(y);
    int indice = index_calc(im.w, im.h, xi, yi, c);
    return im.data[indice];
}

struct lin_transf_coefs{
  float a;
  float b;
};

struct lin_transf_coefs calculate_coefs(int pos1,int pos2){
  float p1 = pos1;
  float p2 = pos2;
  float a = p2 / p1;
  float b = p2 / (2 * p1) - (1.0 / 2.0);
  struct lin_transf_coefs coefs = {a, b};
  return coefs;
}

float lin_transf(struct lin_transf_coefs coef, float x){
  return coef.a * x + coef.b;
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
    image new_im = make_image(w, h, im.c);
    float interp = 0;
    float x, y;
    struct lin_transf_coefs wcoef= calculate_coefs(w, im.w);
    struct lin_transf_coefs hcoef= calculate_coefs(h, im.h);
    for (int i = 0; i < w; i++){
      for (int j = 0; j < h; j++){
          x = lin_transf(wcoef, i);
          y = lin_transf(hcoef, j);
        for (int c = 0; c < im.c; c++){
          interp = nn_interpolate(im, x, y, c);
          new_im.data[c * (w * h) + j * w + i] = interp;
        }
      }
    }
    return new_im;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    int i_floor = floor(x);
    float dx = x - i_floor;
    float val = 0;
    int iprev, ipost, jprev, jpost;
    if (i_floor == -1){
      iprev = 0;
    } else  {
      iprev = i_floor;
    }
    if (i_floor == im.w - 1){
      ipost = im.w - 1;
    } else {
      ipost = i_floor + 1;
    }
    int j_floor = floor(y);
    float dy = y - j_floor;
    if (j_floor == -1){
      jprev = 0;
    } else  {
      jprev = j_floor;
    }
    if (j_floor == im.h - 1){
      jpost = im.h - 1;
    } else {
      jpost = j_floor + 1;
    }
    val += (1 - dx) * (1 - dy) * im.data[index_calc(im.w, im.h, iprev, jprev, c)];
    val += (dx) * (1 - dy) * im.data[index_calc(im.w, im.h, ipost, jprev, c)];
    val += (1 - dx) * (dy) * im.data[index_calc(im.w, im.h, iprev, jpost, c)];
    val += (dx) * (dy) * im.data[index_calc(im.w, im.h, ipost, jpost, c)];
    return val;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    image new_im = make_image(w, h, im.c);
    float interp = 0;
    float x, y;
    struct lin_transf_coefs wcoef= calculate_coefs(w, im.w);
    struct lin_transf_coefs hcoef= calculate_coefs(h, im.h);
    for (int i = 0; i < w; i++){
      for (int j = 0; j < h; j++){
          x = lin_transf(wcoef, i);
          y = lin_transf(hcoef, j);
        for (int c = 0; c < im.c; c++){
          interp = bilinear_interpolate(im, x, y, c);
          new_im.data[c * (w * h) + j * w + i] = interp;
        }
      }
    }
    return new_im;
}

