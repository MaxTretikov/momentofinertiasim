//
//  voxelize.c
//  voxelraster
//
//  Created by Ash Saulsbury on 10/4/20.
//  Copyright © 2020 Ash Saulsbury. All rights reserved.
//


#include "voxellib.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAXBUF 1024
char buffer[MAXBUF];

int main(int argc, const char * argv[]) {

    if (argc != 2) {
      printf("Usage: voxelize <resolution> \n");
      return 0;
    }

    FILE* fp;

    double minf[3], maxf[3];
    double tri1[3], tri2[3], tri3[3];

    Vector3Df_t min, max;

    Vector3Di_t resolution;

    fp = fopen("triangles.txt", "r");

    int line = 0;

    char *p = fgets(buffer, MAXBUF, fp);

    if (sscanf(p, "%lf, %lf, %lf, %lf, %lf, %lf, ",
             &minf[0], &minf[1], &minf[2],
             &maxf[0], &maxf[1], &maxf[2]) != 6);

    set3Df(&max, maxf[0], maxf[1], maxf[2]);
    set3Df(&min, minf[0], minf[1], minf[2]);

    int resi = atoi(argv[1]);

    printf("resolution: %d \n", resi);

    set3Di(&resolution, resi, resi, resi);

    // p = fgets(buffer, MAXBUF, fp);
    //
    // if (sscanf(p, "%lf, %lf, %lf, ",
    //          &resf[0], &resf[1], &resf[3]) != 3);

    init_voxel_box(min, max, resolution);

    do {
      line++;
      p = fgets(buffer, MAXBUF, fp);

      if (p == NULL) break;

      if (sscanf(p, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, ",
               &tri1[0], &tri1[1], &tri1[2],
               &tri2[0], &tri2[1], &tri2[2],
               &tri3[0], &tri3[1], &tri3[2]) != 9) break;

       Triangle_t t;
       set3Df(&(t.vertex[0]), tri1[0], tri1[1], tri1[2]);
       set3Df(&(t.vertex[1]), tri2[0], tri2[1], tri2[2]);
       set3Df(&(t.vertex[2]), tri3[0], tri3[1], tri3[2]);
       voxelizeTriangle(&t);
      } while(1);

      // printf("lines: %d \n", line);

    // for (int z=0; z<resolution.z; z++) {
    //     for (int y=0; y<resolution.y; y++) {
    //         for (int x=0; x<resolution.x; x++) {
    //             int val = getVoxel(x,y,z);
    //             printf("%c",val>0 ? '*':'.');
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    //     // break;
    // }

    FILE *fop;
    fop = fopen("voxels.txt", "a+");

    for (int z=0; z<resolution.z; z++) {
        for (int y=0; y<resolution.y; y++) {
            for (int x=0; x<resolution.x; x++) {
                int val = getVoxel(x,y,z);
                if (val != 0) {
                  fprintf(fop, "%d %d %d\n", x, y, z);
                }
            }
        }
        // break;
    }

    return 0;
}
