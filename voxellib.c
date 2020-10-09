//
//  voxellib.c
//  voxelraster
//
//  Created by Ash Saulsbury on 10/4/20.
//  Copyright © 2020 Ash Saulsbury. All rights reserved.
//

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "voxellib.h"


Vector3Df_t   m_vmin, m_vmax;   // max and min vector ranges
Vector3Di_t   m_vres;      // resolution (integer steps inside voxel box)
int numVoxels;
VoxValue_t *voxelp;

void init_voxel_box(Vector3Df_t min, Vector3Df_t max, Vector3Di_t resolution)
{
    m_vmin = min;
    m_vmax = max;
    m_vres = resolution;

    numVoxels = m_vres.x*m_vres.y*m_vres.z;
    voxelp = calloc(numVoxels, sizeof(voxelp[0]));
    if (voxelp == NULL) {
        fprintf(stderr, "Failed to allocate voxeln");
        exit(1);
    }
}


void clearVoxels()
{
    memset(voxelp, 0, numVoxels*sizeof(voxelp[0]));
}




void set3Df(Vector3Df_t *vecp, float x, float y, float z)
{
    vecp->x = x;
    vecp->y = y;
    vecp->z = z;
}

void set3Di(Vector3Di_t *vecp, int x, int y, int z)
{
    vecp->x = x;
    vecp->y = y;
    vecp->z = z;
}

Vector3Df_t fabs3(Vector3Df_t *fp)
{
    Vector3Df_t temp;
    set3Df(&temp, fabs(fp->x), fabs(fp->y), fabs(fp->z));
    return (temp);
}


// a - b
Vector3Df_t sub3(Vector3Df_t *ap, Vector3Df_t *bp)
{
    Vector3Df_t temp;
    set3Df(&temp, ap->x - bp->x, ap->y - bp->y, ap->z - bp->z);
    return (temp);
}

// a + b
Vector3Df_t add3(Vector3Df_t *ap, Vector3Df_t *bp)
{
    Vector3Df_t temp;
    set3Df(&temp, ap->x + bp->x, ap->y + bp->y, ap->z + bp->z);
    return (temp);
}

// a * scale
Vector3Df_t scale3(Vector3Df_t *ap, float scale)
{
    Vector3Df_t temp;
    set3Df(&temp, ap->x * scale, ap->y * scale, ap->z * scale);
    return (temp);
}


Vector3Df_t cross3(Vector3Df_t *ap, Vector3Df_t *bp)
{
    Vector3Df_t c;
    c.x = (ap->y * bp->z) - (ap->z * bp->y);
    c.y = (ap->z * bp->x) - (ap->x * bp->z);
    c.z = (ap->x * bp->y) - (ap->y * bp->x);
    return (c);
}


float sqlength3(Vector3Df_t *ap)
{
    return (ap->x * ap->x) + (ap->y * ap->y) + (ap->z * ap->z);
}

float length3(Vector3Df_t *ap)
{
    return sqrt(sqlength3(ap));
}

Vector3Df_t normalize3(Vector3Df_t *ap)
{
    Vector3Df_t temp;
    float len = length3(ap);
    set3Df(&temp, ap->x/len, ap->y/len, ap->z/len);
    return (temp);
}


// inverse of a vector.
// it's weird - depends on the definition of vector multiplication,
// but in this case we assume that v/v == unity vector
// which means divide all the elements of v by |v|^2 ...

Vector3Df_t inverse3(Vector3Df_t *ap)
{
    Vector3Df_t temp;
    float val = 1.0 / sqlength3(ap);
    set3Df(&temp, ap->x *val, ap->y * val, ap->z * val);
    return (temp);
}


Vector3Df_t mult3(Vector3Df_t *ap, Vector3Df_t *bp)
{
    Vector3Df_t temp;
    set3Df(&temp, ap->x * bp->x, ap->y * bp->y, ap->z * bp->z);
    return (temp);
}



void setVoxel(int x, int y, int z, uint8_t val)
{
    if (x < 0 || y < 0 || z < 0 || x >= m_vres.x || y >= m_vres.y || z >= m_vres.z) return;
    voxelp[(z*m_vres.y + y)*m_vres.x + x] = val;
}

int getVoxel(int x, int y, int z)
{
    if (x < 0 || y < 0 || z < 0 || x >= m_vres.x || y >= m_vres.y || z >= m_vres.z) return (-1);
    return (voxelp[(z*m_vres.y + y)*m_vres.x + x]);
}




void voxelizeTriangle(Triangle_t *trip)
{
    Vector3Df_t tri[3];

    // transform to unit voxel space
    Vector3Df_t range = sub3(&m_vmax, &m_vmin);

    //    v[i] = (tri[i] - m_vmin) * m_vres / (m_vmax - m_vmin);

    // normalize triangle vertex coords to voxelbox coords
    for (int i=0; i<3; i++) {
        Vector3Df_t vert = sub3(&(trip->vertex[i]), &m_vmin);
        vert.x *= m_vres.x / range.x;
        vert.y *= m_vres.y / range.y;
        vert.z *= m_vres.z / range.z;
        tri[i] = vert;
    }

        // manage this with barycentric coordinates
        // dominate with the longest three sides.

    int p0 = 0;
    int p1 = 1;
    int p2 = 2;

    // printf("P0=%d ; %lf,%lf,%lf\n", p0, tri[p0].x, tri[p0].y, tri[p0].z);
    // printf("P1=%d ; %lf,%lf,%lf\n", p1, tri[p1].x, tri[p1].y, tri[p1].z);
    // printf("P2=%d ; %lf,%lf,%lf\n", p2, tri[p2].x, tri[p2].y, tri[p2].z);

    Vector3Df_t v01 = sub3(&(tri[p1]), &(tri[p0]));
    Vector3Df_t v02 = sub3(&(tri[p2]), &(tri[p0]));
    Vector3Df_t v12 = sub3(&(tri[p2]), &(tri[p1]));
    //
    // printf("v01: %lf,%lf,%lf\n", v01.x, v01.y, v01.z);
    // printf("v02: %lf,%lf,%lf\n", v02.x, v02.y, v02.z);
    // printf("v12: %lf,%lf,%lf\n", v12.x, v12.y, v12.z);

    double len01 = length3(&v01);
    double len02 = length3(&v02);
    double len12 = length3(&v12);

    if (len01 == 0 || len02 == 0 || len12 == 0) return; // line not a triangle

    // where:
    //  x=p2+u.(p0−p2) +v.(p1−p2)
    // rearranges to:
    //  x=u.p0+v.p1+ (1−u−v).p2
    //
    // Following the above we want to walk along the edges p0-p2 and p1-p2 starting from p2
    // for that to work we're hopeing that p0-p2 and p1-p2 are the longest two edges so that
    // p1-p0 is the shortest edge and we don't get gaps as we compute the barycentric coordinates.
    //

    // for our formula p2 is the common vertex, and p0 and p1 are the other two forming edges

    if (len01 <= len02 && len01 <= len12) {
            // v01 and v02 are the two long edges
        p2 = 0; p0 = 2; p1 = 1;
    } else {
            // v01 is not the shortest edge - either v02 or v12 are
        if (len02 <= len12) {
                // v02 is the shortest edge, so v01 and v12 are the two long edges
            p2 = 1; p0 = 0; p1 = 2;
        } else {
                // v12 is the shortest edge, so v01 and v02 are the two long edges
            p2 = 0; p0 = 1; p1 = 2;
        }
    }

    Vector3Df_t vu = sub3(&(tri[p0]), &(tri[p2]));
    Vector3Df_t vv = sub3(&(tri[p1]), &(tri[p2]));
    double stepu = 1.0 / length3(&vu);
    double stepv = 1.0 / length3(&vv);

//    Vector3Df_t roundup;
//    set3Df(&roundup, 0.1, 0.1, 0.1);
    for (double u=0.0; u<1.0; u+=stepu) {
        Vector3Df_t p0p = scale3(&(tri[p0]), u);
        for (double v=0; v<(1.0-u); v+=stepv) {
            double w = 1.0 -u -v;

            Vector3Df_t p1p = scale3(&(tri[p1]), v);
            Vector3Df_t p2p = scale3(&(tri[p2]), w);
            Vector3Df_t p = add3(&p0p, &p1p);
            p = add3(&p, &p2p);
//            p = add3(&p, &roundup);

            setVoxel(p.x, p.y, p.z, 1);
        }
    }
}




//
//
// ----------------------------------------------------------------
// Old code
//

#if 0
void voxelizeTriangle(Triangle_t *trip)
{
    //    clearVoxels();

    float t, tend, vend;
    Vector3Df_t a, b, d1, d2;                // edge variables
    Vector3Df_t c, d, dstep, mask, s;        // 2DDA variables
    Vector3Df_t v[3];

    // transform to unit voxel space
    Vector3Df_t range = sub3(&m_vmax, &m_vmin);

    //    v[i] = (tri[i] - m_vmin) * m_vres / (m_vmax - m_vmin);

    // normalize triangle vertex coords to voxelbox coords
    for (int i=0; i<3; i++) {
        Vector3Df_t vert = sub3(&(trip->vertex[i]), &m_vmin);
        vert.x *= m_vres.x / range.x;
        vert.y *= m_vres.y / range.y;
        vert.z *= m_vres.z / range.z;
        v[i] = vert;
    }

    //    d = d.Cross(v[1] - v[0], v[2] - v[0]);
    Vector3Df_t va = sub3(&(v[1]), &(v[0]));
    Vector3Df_t vb = sub3(&(v[2]), &(v[0]));
    d = cross3(&va, &vb);
    d = normalize3(&d);

    // using the surface normal of the triangle to determine which axis it mostly faces towards

    // determine raster axis
    int axis = fabs(d.x) > fabs(d.y) ? (fabs(d.x) > fabs(d.z) ? 0 : 2) : (fabs(d.y) > fabs(d.z)) ? 1 : 2;

    // sort vertices
    int p0 = 0;
    int p1 = 1;
    int p2 = 2;
    int p3; // temp
    if (axis == 1) {
        // use x axis (triangle oriented toward y)
        if (v[p0].x > v[p1].x) { p3 = p0; p0 = p1; p1 = p3; }
        if (v[p0].x > v[p2].x) { p3 = p0; p0 = p2; p2 = p3; }
        if (v[p1].x > v[p2].x) { p3 = p1; p1 = p2; p2 = p3; }

        // edge setup
        a = v[p0];
        b = v[p0];
        //        d1 = (v[p1] - v[p0]) / (v[p1].x - v[p0].x);
        //        d2 = (v[p2] - v[p0]) / (v[p2].x - v[p0].x);
        d1 = sub3(&(v[p1]), &(v[p0]));
        d1 = scale3(&d1, 1.0 / (v[p1].x - v[p0].x));
        d2 = sub3(&(v[p2]), &(v[p0]));
        d2 = scale3(&d2, 1.0 / (v[p2].x - v[p0].x));


        //        a = a + d1 * float(1.0 - v[p0].x + floor(v[p0].x));
        //        b = b + d2 * float(1.0 - v[p0].x + floor(v[p0].x));
        d1 = scale3(&d1, 1.0 - v[p0].x + floor(v[p0].x));
        a = add3(&a, &d1);
        d2 = scale3(&d2, 1.0 - v[p0].x + floor(v[p0].x));
        b = add3(&b, &d2);
        vend = v[p1].x;

        // rasterize - outer loop is move along x axis
        for (; b.x < v[p2].x; ) {

            // prepare 2DDA
            t = 0;
            c = a;
            d = sub3(&b, &c);
            tend = length3(&d);
            d = normalize3(&d);        // scan vector
            d.x = 1;
            set3Df(&dstep,
                   (d.x > 0) ? 1 : -1,
                   (d.y > 0) ? 1 : -1,
                   (d.z > 0) ? 1 : -1);    // signed direction
            d = fabs3(&d);
            set3Df(&s, 0,
                   ((((int)c.y) - c.y + 0.5f)*dstep.y + 0.5) / d.y + t,
                   ((((int)c.z) - c.z + 0.5f)*dstep.z + 0.5) / d.z + t);

            Vector3Df_t dinv = inverse3(&d);

            // 2DDA in the YZ-plane
            while (t < tend) {
                setVoxel(c.x + 0.1, c.y, c.z, 1);
                setVoxel(c.x - 0.1, c.y, c.z, 1);
                //                mask = (s.y < s.z) ? Vector3DI(0, 1, 0) : Vector3DI(0, 0, 1);    // choose next voxel
                //                t = mask.y ? s.y : s.z;                // advance t
                //                s += mask / d;                        // advance x/y/z intercepts by the inverse normal (not obvious)
                //                c += mask * dstep;                    // advance to next voxel
                if (s.y < s.z) {
                    set3Df(&mask, 0, 1, 0);
                } else {
                    set3Df(&mask, 0, 0, 1);
                } // choose next voxel
                t = mask.y ? s.y : s.z;                // advance t

                //                s += mask / d;                        // advance x/y/z intercepts by the inverse normal (not obvious)
                //                c += mask * dstep;                    // advance to next voxel
                Vector3Df_t step = mult3(&mask, &dinv);
                s = add3(&s, &step); // advance x/y/z intercepts by the inverse normal (not obvious)
                step = mult3(&mask, &dstep);        // re-use step variable
                c = add3(&mask, &dstep);                    // advance to next voxel
            }
            c = b;
            setVoxel(c.x + 0.1, c.y, c.z, 1);
            setVoxel(c.x - 0.1, c.y, c.z, 1);

            // advance along edges
            a = add3(&a, &d1);
            b = add3(&b, &d2);
            if (a.x >= vend) {            // transition for third edge
                vend = a.x - v[p1].x;    // measure amount we overshot

                Vector3Df_t sub = scale3(&d1, vend);
                a = sub3(&a, &sub);            // back up to corner

                //                d1 = (v[p2] - v[p1]) / (v[p2].x - v[p1].x);
                //                a += d1*vend;            // advance to y-unit along new edge
                Vector3Df_t add = sub3(&(v[p2]), &(v[p1]));
                add = scale3(&add, vend / (v[p2].x - v[p1].x));
                a = add3(&a, &add);            // advance to y-unit along new edge
                vend = v[p2].x;
            }
        }


    } else {

        printf("P0=%d ; %lf,%lf,%lf\n", p0, v[p0].x, v[p0].y, v[p0].z);
        printf("P1=%d ; %lf,%lf,%lf\n", p1, v[p1].x, v[p1].y, v[p1].z);
        printf("P2=%d ; %lf,%lf,%lf\n", p2, v[p2].x, v[p2].y, v[p2].z);

        // use y axis (triangle oriented toward x or z)
        if (v[p0].y > v[p1].y) { p3 = p0; p0 = p1; p1 = p3; }
        if (v[p0].y > v[p2].y) { p3 = p0; p0 = p2; p2 = p3; }
        if (v[p1].y > v[p2].y) { p3 = p1; p1 = p2; p2 = p3; }

        printf("\n");
        printf("P0=%d ; %lf,%lf,%lf\n", p0, v[p0].x, v[p0].y, v[p0].z);
        printf("P1=%d ; %lf,%lf,%lf\n", p1, v[p1].x, v[p1].y, v[p1].z);
        printf("P2=%d ; %lf,%lf,%lf\n", p2, v[p2].x, v[p2].y, v[p2].z);


        // edge setup
        //
        a = v[p0];
        b = v[p0];


        //        d1 = (v[p1] - v[p0]) / (v[p1].y - v[p0].y);
        //        d2 = (v[p2] - v[p0]) / (v[p2].y - v[p0].y);
        d1 = sub3(&(v[p1]), &(v[p0]));
        d1 = scale3(&d1, 1.0 / (v[p1].y - v[p0].y));
        d2 = sub3(&(v[p2]), &(v[p0]));
        d2 = scale3(&d2, 1.0 / (v[p2].y - v[p0].y));

        //        a = a + d1 * float(1.0 - v[p0].y + floor(v[p0].y));
        //        b = b + d2 * float(1.0 - v[p0].y + floor(v[p0].y));
        d1 = scale3(&d1, (float)(1.0 - v[p0].y + floor(v[p0].y)));
        a = add3(&a, &d1);
        d2 = scale3(&d2, (float)(1.0 - v[p0].y + floor(v[p0].y)));
        b = add3(&b, &d2);
        vend = v[p1].y;


        // rasterize along y axis
        for (; b.y < v[p2].y; ) {

            // prepare 2DDA
            t = 0;
            c = a;
            d = sub3(&b, &c);
            tend = length3(&d);
            d = normalize3(&d);        // scan vector
            d.y = 1;
            set3Df(&dstep,
                   (d.x > 0) ? 1 : -1,
                   (d.y > 0) ? 1 : -1,
                   (d.z > 0) ? 1 : -1);    // signed direction
            d = fabs3(&d);
            set3Df(&s,
                   (((int)(c.x) - c.x + 0.5f)*dstep.x + 0.5) / d.x + t,
                   0,
                   (((int)(c.z) - c.z + 0.5f)*dstep.z + 0.5) / d.z + t);

            Vector3Df_t dinv = inverse3(&d);

            // 2DDA in the XZ-plane
            while (t < tend) {
                setVoxel(c.x, c.y + 0.1, c.z, 1);
                setVoxel(c.x, c.y - 0.1, c.z, 1);

                //                mask = (s.x < s.z) ? Vector3DI(1, 0, 0) : Vector3DI(0, 0, 1);    // choose next voxel
                if (s.x < s.z) {
                    set3Df(&mask, 1, 0, 0);
                } else {
                    set3Df(&mask, 0, 0, 1);
                }
                t = mask.x ? s.x : s.z;                // advance t
                //                s += mask / d;                        // advance x/y/z intercepts by the inverse normal (not obvious)
                //                c += mask * dstep;                    // advance to next voxel
                Vector3Df_t step = mult3(&mask, &dinv);
                c = add3(&c, &step);
            }
            c = b;
            setVoxel(c.x, c.y + 0.1, c.z, 1);
            setVoxel(c.x, c.y - 0.1, c.z, 1);

            // advance along edges
            a = add3(&a, &d1);
            b = add3(&b, &d2);
            if (a.y >= vend) {            // transition for third edge
                vend = a.y - v[p1].y;    // measure amount we overshot
                //                a -= d1*vend;            // back up to corner
                Vector3Df_t add = scale3(&d1, vend);
                a = sub3(&a, &add);

                add = sub3(&(v[p2]), &(v[p1]));
                //                d1 = (v[p2] - v[p1]) / (v[p2].y - v[p1].y);
                //                a += d1*vend;            // advance to y-unit along new edge
                add = scale3(&add, vend / (v[p2].y - v[p1].y));
                a = add3(&a, &add);            // advance to y-unit along new edge
                vend = v[p2].y;
            }
        }
    }

    //    m_voxel_visit = m_voxel_set;
}

#endif
