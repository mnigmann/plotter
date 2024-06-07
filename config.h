#ifndef CONFIG_H
#define CONFIG_H

#define INIT_WIDTH 500
#define INIT_HEIGHT 500
#define INIT_SCALE 25
#define COLOR_AXES 64
#define COLOR_MAJOR_GRID 160
#define COLOR_MINOR_GRID 220
#define COLOR_SPECIAL_POINT 0xdcdcdc

// Distance between axis and number, in pixels
#define AXIS_DISTANCE 3
// Maximum distance between axis and cursor when scaling axes independently
#define AXIS_MARGIN 20

#define POINT_SIZE 3
#define MAX_POLYGON_SIZE 4

// in scale units
#define PLOT_DT_DERIV 1e-7
// in pixels
#define PLOT_MAX_ARC_LENGTH 3
#define PLOT_MIN_DT 1e-3
#define PLOT_IMPLICIT_MAXDEPTH 6
#define PLOT_IMPLICIT_HARDMAX 8
#define PLOT_USE_INEQUALITY
#define PLOT_EXTREMUM_MAXITER 10
#define PLOT_EXTREMUM_EPS 1e-8
#define PLOT_DISCONTINUITY_THRESHOLD 4
#define PLOT_CONTINUITY_THRESHOLD 0.8
#define PLOT_DISCONTINUITY_MAXITER 16
#define CLICK_MAX_RADIUS 10
#define CLICK_MAX_ITER 4
#define CLICK_IMPLICIT_MAX_ERROR 1e-6

#define INTEGRATE_DX 1e-2

#define PI 3.14159265359
#define SCROLL_AMOUNT 1.1
#define TICK_SIZE 125

#define CACHE_BLOCK_SIZE 256

#define DEBUG_EVAL
#define DEBUG_PLOT

#endif
