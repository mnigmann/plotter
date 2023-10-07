#include <stdint.h>
#include <stdio.h>
#include <string.h>


enum { XY_SHIFT = 16, XY_ONE = 1 << XY_SHIFT, DRAWING_STORAGE_BLOCK = (1<<12) - 256 };


/* Correction table depent on the slope */
const uint8_t SlopeCorrTable[] = {
    181, 181, 181, 182, 182, 183, 184, 185, 187, 188, 190, 192, 194, 196, 198, 201,
    203, 206, 209, 211, 214, 218, 221, 224, 227, 231, 235, 238, 242, 246, 250, 254
};

/* Gaussian for antialiasing filter */
const uint8_t FilterTable[] = {
    168, 177, 185, 194, 202, 210, 218, 224, 231, 236, 241, 246, 249, 252, 254, 254,
    254, 254, 252, 249, 246, 241, 236, 231, 224, 218, 210, 202, 194, 185, 177, 168,
    158, 149, 140, 131, 122, 114, 105, 97, 89, 82, 75, 68, 62, 56, 50, 45,
    40, 36, 32, 28, 25, 22, 19, 16, 14, 12, 11, 9, 8, 7, 5, 5
};


uint8_t clipLine(uint64_t size_width, uint64_t size_height, int64_t *pt1x, int64_t *pt1y, int64_t *pt2x, int64_t *pt2y) {
    int c1, c2;
    int64_t right = size_width-1, bottom = size_height-1;

    if( size_width <= 0 || size_height <= 0 )
        return 0;

    c1 = (*pt1x < 0) + (*pt1x > right) * 2 + (*pt1y < 0) * 4 + (*pt1y > bottom) * 8;
    c2 = (*pt2x < 0) + (*pt2x > right) * 2 + (*pt2y < 0) * 4 + (*pt2y > bottom) * 8;

    if ((c1 & c2) == 0 && (c1 | c2) != 0) {
        int64_t a;
        if (c1 & 12) {
            a = c1 < 8 ? 0 : bottom;
            *pt1x += (int64_t)((double)(a - *pt1y) * (*pt2x - *pt1x) / (*pt2y - *pt1y));
            *pt1y = a;
            c1 = (*pt1x < 0) + (*pt1x > right) * 2;
        }
        if (c2 & 12) {
            a = c2 < 8 ? 0 : bottom;
            *pt2x += (int64_t)((double)(a - *pt2y) * (*pt2x - *pt1x) / (*pt2y - *pt1y));
            *pt2y = a;
            c2 = (*pt2x < 0) + (*pt2x > right) * 2;
        }
        if ((c1 & c2) == 0 && (c1 | c2) != 0) {
            if (c1) {
                a = c1 == 1 ? 0 : right;
                *pt1y += (int64_t)((double)(a - *pt1x) * (*pt2y - *pt1y) / (*pt2x - *pt1x));
                *pt1x = a;
                c1 = 0;
            }
            if (c2) {
                a = c2 == 1 ? 0 : right;
                *pt2y += (int64_t)((double)(a - *pt2x) * (*pt2y - *pt1y) / (*pt2x - *pt1x));
                *pt2x = a;
                c2 = 0;
            }
        }
    }

    return (c1 | c2) == 0;
}


void LineAA(uint8_t *ptr, uint64_t size_width, uint64_t size_height, int64_t pt1x, int64_t pt1y, int64_t pt2x, int64_t pt2y, const void* color ) {
    int64_t dx, dy;
    int ecount, scount = 0;
    int slope;
    int64_t ax, ay;
    int64_t x_step, y_step;
    int64_t i, j;
    int ep_table[9];
    int cb = ((uint8_t*)color)[0], cg = ((uint8_t*)color)[1], cr = ((uint8_t*)color)[2], ca = ((uint8_t*)color)[3];
    int _cb, _cg, _cr, _ca;
    int nch = 3;
    uint64_t step = nch * size_width;
    //pt1x <<= XY_SHIFT; pt1y <<= XY_SHIFT; pt2x <<= XY_SHIFT; pt2y <<= XY_SHIFT;

    /*if( !((nch == 1 || nch == 3 || nch == 4) && img.depth() == CV_8U) )
    {
        Line(img, Point((int)(pt1x>>XY_SHIFT), (int)(pt1y>>XY_SHIFT)), Point((int)(pt2x>>XY_SHIFT), (int)(pt2y>>XY_SHIFT)), color);
        return;
    }*/

    //size_width <<= XY_SHIFT;
    //size_height <<= XY_SHIFT;
    //    return;
    int64_t large_width = size_width << XY_SHIFT;
    int64_t large_height = size_height << XY_SHIFT;
    //if ((pt1x < 0) || (pt1y < 0) || (pt2x < 0) || (pt2y < 0) || (pt1x >= large_width) || (pt2y >= large_height) || (pt2x >= large_width) || (pt2y >= large_height))
    //if (((pt1x < 0) && (pt2x < 0)) || ((pt1y < 0) && (pt2y < 0)) || ((pt1x >= large_width) && (pt2x >= large_width)) || ((pt1y >= large_height) && (pt2y >= large_height)))
    if (!clipLine(large_width, large_height, &pt1x, &pt1y, &pt2x, &pt2y))
        return;

    dx = pt2x - pt1x;
    dy = pt2y - pt1y;

    j = dx < 0 ? -1 : 0;
    ax = (dx ^ j) - j;      // abs(dx)
    i = dy < 0 ? -1 : 0;
    ay = (dy ^ i) - i;      // abs(dy)

    if( ax > ay )
    {
        dy = (dy ^ j) - j;  // negate dy if dx < 0
        pt1x ^= pt2x & j;   // swap pt1 and pt2 if dx < 0
        pt2x ^= pt1x & j;
        pt1x ^= pt2x & j;
        pt1y ^= pt2y & j;
        pt2y ^= pt1y & j;
        pt1y ^= pt2y & j;

        x_step = XY_ONE;
        y_step = (int64_t)((uint64_t)dy << XY_SHIFT) / (ax | 1);
        pt2x += XY_ONE;
        ecount = (int)((pt2x >> XY_SHIFT) - (pt1x >> XY_SHIFT));
        j = -(pt1x & (XY_ONE - 1));
        pt1y += ((y_step * j) >> XY_SHIFT) + (XY_ONE >> 1);
        slope = (y_step >> (XY_SHIFT - 5)) & 0x3f;
        slope ^= (y_step < 0 ? 0x3f : 0);

        /* Get 4-bit fractions for end-point adjustments */
        i = (pt1x >> (XY_SHIFT - 7)) & 0x78;
        j = (pt2x >> (XY_SHIFT - 7)) & 0x78;
    }
    else
    {
        dx = (dx ^ i) - i;
        pt1x ^= pt2x & i;
        pt2x ^= pt1x & i;
        pt1x ^= pt2x & i;
        pt1y ^= pt2y & i;
        pt2y ^= pt1y & i;
        pt1y ^= pt2y & i;

        x_step = (int64_t)((uint64_t)dx << XY_SHIFT) / (ay | 1);
        y_step = XY_ONE;
        pt2y += XY_ONE;
        ecount = (int)((pt2y >> XY_SHIFT) - (pt1y >> XY_SHIFT));
        j = -(pt1y & (XY_ONE - 1));
        pt1x += ((x_step * j) >> XY_SHIFT) + (XY_ONE >> 1);
        slope = (x_step >> (XY_SHIFT - 5)) & 0x3f;
        slope ^= (x_step < 0 ? 0x3f : 0);

        /* Get 4-bit fractions for end-point adjustments */
        i = (pt1y >> (XY_SHIFT - 7)) & 0x78;
        j = (pt2y >> (XY_SHIFT - 7)) & 0x78;
    }

    slope = (slope & 0x20) ? 0x100 : SlopeCorrTable[slope];

    /* Calc end point correction table */
    {
        int t0 = slope << 7;
        int t1 = ((0x78 - (int)i) | 4) * slope;
        int t2 = ((int)j | 4) * slope;

        ep_table[0] = 0;
        ep_table[8] = slope;
        ep_table[1] = ep_table[3] = ((((j - i) & 0x78) | 4) * slope >> 8) & 0x1ff;
        ep_table[2] = (t1 >> 8) & 0x1ff;
        ep_table[4] = ((((j - i) + 0x80) | 4) * slope >> 8) & 0x1ff;
        ep_table[5] = ((t1 + t0) >> 8) & 0x1ff;
        ep_table[6] = (t2 >> 8) & 0x1ff;
        ep_table[7] = ((t2 + t0) >> 8) & 0x1ff;
    }

    if( nch == 3 )
    {
        #define  ICV_PUT_POINT(x, y)        \
        {                                   \
            uint8_t* tptr = ptr + (x)*3 + (y)*step; \
            _cb = tptr[0];                  \
            _cb += ((cb - _cb)*a + 127)>> 8;\
            _cb += ((cb - _cb)*a + 127)>> 8;\
            _cg = tptr[1];                  \
            _cg += ((cg - _cg)*a + 127)>> 8;\
            _cg += ((cg - _cg)*a + 127)>> 8;\
            _cr = tptr[2];                  \
            _cr += ((cr - _cr)*a + 127)>> 8;\
            _cr += ((cr - _cr)*a + 127)>> 8;\
            tptr[0] = (uint8_t)_cb;           \
            tptr[1] = (uint8_t)_cg;           \
            tptr[2] = (uint8_t)_cr;           \
        }
        if( ax > ay )       // Line is wider than it is tall
        {
            int x = (int)(pt1x >> XY_SHIFT);

            for( ; ecount >= 0; x++, pt1y += y_step, scount++, ecount-- )
            {
                if (((unsigned)x >= (unsigned)size_width) || (x < 0))
                    continue;
                int y = (int)((pt1y >> XY_SHIFT) - 1);

                int ep_corr = ep_table[(((scount >= 2) + 1) & (scount | 2)) * 3 +
                                       (((ecount >= 2) + 1) & (ecount | 2))];
                int a, dist = (pt1y >> (XY_SHIFT - 5)) & 31;
#ifdef DEBUG
                printf("ecount: %d, scount: %d, ep_corr: %d, dist: %d, pt1y: %08x\n", ecount, scount, ep_corr, dist, pt1y);
#endif

                a = (ep_corr * FilterTable[dist + 32] >> 8) & 0xff;
                if( (unsigned)y < (unsigned)size_height )
                    ICV_PUT_POINT(x, y)

                a = (ep_corr * FilterTable[dist] >> 8) & 0xff;
                if( (unsigned)(y+1) < (unsigned)size_height )
                    ICV_PUT_POINT(x, y+1)

                a = (ep_corr * FilterTable[63 - dist] >> 8) & 0xff;
                if( (unsigned)(y+2) < (unsigned)size_height )
                    ICV_PUT_POINT(x, y+2)
            }
        }
        else
        {
            int y = (int)(pt1y >> XY_SHIFT);

            for( ; ecount >= 0; y++, pt1x += x_step, scount++, ecount-- )
            {
                if (((unsigned)y >= (unsigned)size_height) || (y < 0))
                    continue;
                int x = (int)((pt1x >> XY_SHIFT) - 1);
                int ep_corr = ep_table[(((scount >= 2) + 1) & (scount | 2)) * 3 +
                                       (((ecount >= 2) + 1) & (ecount | 2))];
                int a, dist = (pt1x >> (XY_SHIFT - 5)) & 31;

                a = (ep_corr * FilterTable[dist + 32] >> 8) & 0xff;
                if( (unsigned)x < (unsigned)size_width )
                    ICV_PUT_POINT(x, y)

                a = (ep_corr * FilterTable[dist] >> 8) & 0xff;
                if( (unsigned)(x+1) < (unsigned)size_width )
                    ICV_PUT_POINT(x+1, y)

                a = (ep_corr * FilterTable[63 - dist] >> 8) & 0xff;
                if( (unsigned)(x+2) < (unsigned)size_width )
                    ICV_PUT_POINT(x+2, y)
            }
        }
        #undef ICV_PUT_POINT
    }
    else if(nch == 1)
    {
        #define ICV_PUT_POINT(x, y)         \
        {                                   \
            uint8_t* tptr = ptr + (x) + (y) * step; \
            _cb = tptr[0];                  \
            _cb += ((cb - _cb)*a + 127)>> 8;\
            _cb += ((cb - _cb)*a + 127)>> 8;\
            tptr[0] = (uint8_t)_cb;           \
        }

        if( ax > ay )
        {
            int x = (int)(pt1x >> XY_SHIFT);

            for( ; ecount >= 0; x++, pt1y += y_step, scount++, ecount-- )
            {
                if( (unsigned)x >= (unsigned)size_width )
                    continue;
                int y = (int)((pt1y >> XY_SHIFT) - 1);

                int ep_corr = ep_table[(((scount >= 2) + 1) & (scount | 2)) * 3 +
                                       (((ecount >= 2) + 1) & (ecount | 2))];
                int a, dist = (pt1y >> (XY_SHIFT - 5)) & 31;

                a = (ep_corr * FilterTable[dist + 32] >> 8) & 0xff;
                if( (unsigned)y < (unsigned)size_height )
                    ICV_PUT_POINT(x, y)

                a = (ep_corr * FilterTable[dist] >> 8) & 0xff;
                if( (unsigned)(y+1) < (unsigned)size_height )
                    ICV_PUT_POINT(x, y+1)

                a = (ep_corr * FilterTable[63 - dist] >> 8) & 0xff;
                if( (unsigned)(y+2) < (unsigned)size_height )
                    ICV_PUT_POINT(x, y+2)
            }
        }
        else
        {
            int y = (int)(pt1y >> XY_SHIFT);

            for( ; ecount >= 0; y++, pt1x += x_step, scount++, ecount-- )
            {
                if( (unsigned)y >= (unsigned)size_height )
                    continue;
                int x = (int)((pt1x >> XY_SHIFT) - 1);
                int ep_corr = ep_table[(((scount >= 2) + 1) & (scount | 2)) * 3 +
                                       (((ecount >= 2) + 1) & (ecount | 2))];
                int a, dist = (pt1x >> (XY_SHIFT - 5)) & 31;

                a = (ep_corr * FilterTable[dist + 32] >> 8) & 0xff;
                if( (unsigned)x < (unsigned)size_width )
                    ICV_PUT_POINT(x, y)

                a = (ep_corr * FilterTable[dist] >> 8) & 0xff;
                if( (unsigned)(x+1) < (unsigned)size_width )
                    ICV_PUT_POINT(x+1, y)

                a = (ep_corr * FilterTable[63 - dist] >> 8) & 0xff;
                if( (unsigned)(x+2) < (unsigned)size_width )
                    ICV_PUT_POINT(x+2, y)
            }
        }
        #undef ICV_PUT_POINT
    }
    else
    {
        #define  ICV_PUT_POINT(x, y)        \
        {                                   \
            uint8_t* tptr = ptr + (x)*4 + (y)*step; \
            _cb = tptr[0];                  \
            _cb += ((cb - _cb)*a + 127)>> 8;\
            _cb += ((cb - _cb)*a + 127)>> 8;\
            _cg = tptr[1];                  \
            _cg += ((cg - _cg)*a + 127)>> 8;\
            _cg += ((cg - _cg)*a + 127)>> 8;\
            _cr = tptr[2];                  \
            _cr += ((cr - _cr)*a + 127)>> 8;\
            _cr += ((cr - _cr)*a + 127)>> 8;\
            _ca = tptr[3];                  \
            _ca += ((ca - _ca)*a + 127)>> 8;\
            _ca += ((ca - _ca)*a + 127)>> 8;\
            tptr[0] = (uint8_t)_cb;           \
            tptr[1] = (uint8_t)_cg;           \
            tptr[2] = (uint8_t)_cr;           \
            tptr[3] = (uint8_t)_ca;           \
        }
        if( ax > ay )
        {
            int x = (int)(pt1x >> XY_SHIFT);

            for( ; ecount >= 0; x++, pt1y += y_step, scount++, ecount-- )
            {
                if( (unsigned)x >= (unsigned)size_width )
                    continue;
                int y = (int)((pt1y >> XY_SHIFT) - 1);

                int ep_corr = ep_table[(((scount >= 2) + 1) & (scount | 2)) * 3 +
                                       (((ecount >= 2) + 1) & (ecount | 2))];
                int a, dist = (pt1y >> (XY_SHIFT - 5)) & 31;

                a = (ep_corr * FilterTable[dist + 32] >> 8) & 0xff;
                if( (unsigned)y < (unsigned)size_height )
                    ICV_PUT_POINT(x, y)

                a = (ep_corr * FilterTable[dist] >> 8) & 0xff;
                if( (unsigned)(y+1) < (unsigned)size_height )
                    ICV_PUT_POINT(x, y+1)

                a = (ep_corr * FilterTable[63 - dist] >> 8) & 0xff;
                if( (unsigned)(y+2) < (unsigned)size_height )
                    ICV_PUT_POINT(x, y+2)
            }
        }
        else
        {
            int y = (int)(pt1y >> XY_SHIFT);

            for( ; ecount >= 0; y++, pt1x += x_step, scount++, ecount-- )
            {
                if( (unsigned)y >= (unsigned)size_height )
                    continue;
                int x = (int)((pt1x >> XY_SHIFT) - 1);
                int ep_corr = ep_table[(((scount >= 2) + 1) & (scount | 2)) * 3 +
                                       (((ecount >= 2) + 1) & (ecount | 2))];
                int a, dist = (pt1x >> (XY_SHIFT - 5)) & 31;

                a = (ep_corr * FilterTable[dist + 32] >> 8) & 0xff;
                if( (unsigned)x < (unsigned)size_width )
                    ICV_PUT_POINT(x, y)

                a = (ep_corr * FilterTable[dist] >> 8) & 0xff;
                if( (unsigned)(x+1) < (unsigned)size_width )
                    ICV_PUT_POINT(x+1, y)

                a = (ep_corr * FilterTable[63 - dist] >> 8) & 0xff;
                if( (unsigned)(x+2) < (unsigned)size_width )
                    ICV_PUT_POINT(x+2, y)
            }
        }
        #undef ICV_PUT_POINT
    }
}

#ifdef DEBUG
uint8_t data[2700];
uint8_t color[4] = {255, 255, 255, 255};

int main() {
    memset(data, 0, 2700);
    LineAA(data, 30, 30, 4<<16, 4<<16, 20<<16, 12<<16, color);

    for (int i=0; i < 30; i++) {
        for (int j=0; j < 30; j++) {
            printf("%4d", data[3*(30*i+j)]);
        }
        printf("\n");
    }
}
#endif

