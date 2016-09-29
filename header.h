#ifndef TEMA2_SD_HEADER_H
#define TEMA2_SD_HEADER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef struct Pixel {
    unsigned char colors[4];
} Pixel ;
typedef struct QuadtreeNode {
    unsigned char blue, green, red, reserved;
    uint32_t area;
    int32_t top_left, top_right;
    int32_t bottom_left, bottom_right;
} __attribute__ ((packed)) QuadtreeNode;

typedef struct Tree {
    unsigned char blue, green, red, reserved;
    uint32_t area;
    struct Tree *top_left, *top_right, *bottom_left, *bottom_right;
} Tree;

int read(Pixel ***, FILE *, FILE *);
void buildQuadtree(Pixel **, Tree *,  int ,  int ,  int , int );

#endif //TEMA2_SD_HEADER_H
