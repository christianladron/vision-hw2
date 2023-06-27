#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    if (x <0) x = 0;
    if (y <0) y = 0;
    if (x >= im.w) x = im.w - 1;
    if (y >= im.h) y = im.h - 1;
    float v = im.data[c * (im.w * im.h) + y * im.w + x];
    return v;
}

void set_pixel(image im, int x, int y, int c, float v)
{
    if (x <0 || y <0 || c <0 || x >= im.w || y >= im.h || c>= im.c) return;
    im.data[c * (im.w * im.h) + y * im.w + x] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    for (int i = 0; i < im.w; i++){
    for (int j = 0; j < im.h; j++){
    for (int k = 0; k < im.c; k++){
      float v = get_pixel(im, i, j, k);
      set_pixel(copy, i, j, k, v);
    }}}
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    for (int x = 0; x < im.w; x++){
    for (int y = 0; y < im.h; y++){
      float v = 0;
      v += 0.299 * get_pixel(im, x, y, 0);
      v += 0.587 * get_pixel(im, x, y, 1);
      v += 0.114 * get_pixel(im, x, y, 2);
      set_pixel(gray, x, y, 0, v);
    }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    for (int i = 0; i < im.w; i++){
    for (int j = 0; j < im.h; j++){
      float vi = get_pixel(im, i, j, c);
      vi = vi + v;
      set_pixel(im, i, j, c, vi);
    }}
}

void clamp_image(image im)
{
    for (int i = 0; i < im.w; i++){
    for (int j = 0; j < im.h; j++){
    for (int k = 0; k < im.c; k++){
      float v = get_pixel(im, i, j, k);
      if (v <0) set_pixel(im, i, j, k, 0);
      else if (v > 1) set_pixel(im, i, j, k, 1);
    }}}
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    for (int i = 0; i < im.w; i++){
    for (int j = 0; j < im.h; j++){
      float hp, h;
      float r = get_pixel(im, i, j, 0);
      float g = get_pixel(im, i, j, 1);
      float b = get_pixel(im, i, j, 2);
      float v = three_way_max(r, g, b);
      float m = three_way_min(r, g, b);
      float s = 0;
      float C = (v - m);
      if (v > 0) s = C / v;
      if (C == 0) hp = 0;
      else if (v == r){
        hp = (g - b) / C;
      }
      else if (v == g){
        hp = (b - r) / C + 2;
      }
      else if (v == b){
        hp = (r - g) / C + 4;
      }
      if (hp < 0){
        h = hp / 6 + 1;
      }
      else{
        h = hp / 6;
      }
      set_pixel(im, i, j, 0, h);
      set_pixel(im, i, j, 1, s);
      set_pixel(im, i, j, 2, v);
    }}
}

void hsv_to_rgb(image im)
{
    for (int i = 0; i < im.w; i++){
    for (int j = 0; j < im.h; j++){
      float r, g, b;
      float h = get_pixel(im, i, j, 0);
      float s = get_pixel(im, i, j, 1);
      float v = get_pixel(im, i, j, 2);
      float c = v * s;
      float hp = h * 6;
      float x = c * (1 - fabs(fmod(hp, 2.0) - 1));
      if (hp < 1){
        r = c;
        g = x;
        b = 0;
      }
      else if (hp < 2){
        r = x;
        g = c;
        b = 0;
      }
      else if (hp < 3){
        r = 0;
        g = c;
        b = x;
      }
      else if (hp < 4){
        r = 0;
        g = x;
        b = c;
      }
      else if (hp < 5){
        r = x;
        g = 0;
        b = c;
      }
      else if (hp < 6){
        r = c;
        g = 0;
        b = x;
      }
      float m = v - c;
      r += m;
      g += m;
      b += m;
      set_pixel(im, i, j, 0, r);
      set_pixel(im, i, j, 1, g);
      set_pixel(im, i, j, 2, b);
    }}
}

void scale_image(image im, int c, float v)
{
    for (int i = 0; i < im.w; i++){
    for (int j = 0; j < im.h; j++){
      float vi = get_pixel(im, i, j, c);
      vi = vi * v;
      set_pixel(im, i, j, c, vi);
    }}
}

