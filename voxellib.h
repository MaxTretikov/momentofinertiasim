//
//  voxellib.h
//  voxelraster
//
//  Created by Ash Saulsbury on 10/4/20.
//  Copyright © 2020 Ash Saulsbury. All rights reserved.
//

#ifndef voxellib_h
#define voxellib_h

#include <stdio.h>
#include <stdint.h>
#include <math.h>


typedef struct {
    float x,y,z;
} Vector3Df_t;

typedef struct {
    int x,y,z;
} Vector3Di_t;

typedef struct {
    Vector3Df_t vertex[3];
} Triangle_t;

typedef uint8_t VoxValue_t;


extern void init_voxel_box(Vector3Df_t min, Vector3Df_t max, Vector3Di_t resolution);
extern void clearVoxels(void);
extern void set3Df(Vector3Df_t *vecp, float x, float y, float z);
extern void set3Di(Vector3Di_t *vecp, int x, int y, int z);
extern void voxelizeTriangle(Triangle_t *trip);
extern void setVoxel(int x, int y, int z, uint8_t val);
extern int getVoxel(int x, int y, int z);

#endif /* voxellib_h */
