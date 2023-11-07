#include <gtk/gtk.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "parse.h"
#include "functions.h"
#include "linalg_functions.h"
#include "config.h"
#include "treeview.h"

#define SCALE_X(v) ((uint16_t)(250*(1+v)))
#define SCALE_Y(v) ((uint16_t)(250*(1-v)))

#define SCALE_XF(v) ((int64_t)(250L*(1<<16)*(1+v)))
#define SCALE_YF(v) ((int64_t)(250L*(1<<16)*(1-v)))

#define SCALE_XK(v) ((v - window_x0)*xscale)
#define SCALE_YK(v) ((window_y1 - v)*yscale)

#define CLIP_WIDTH(v) (v<0 ? 0 : (v>WIDTH ? WIDTH : v))
#define CLIP_HEIGHT(v) (v<0 ? 0 : (v>HEIGHT ? HEIGHT : v))
#define IN_BOUNDS(x, y) ((window_x0 <= x) && (x <= window_x1) && (window_y0 <= y) && (y <= window_y1))
#define SET_COLOR(cr, color) cairo_set_source_rgb(cr, color[0]/256.0, color[1]/256.0, color[2]/256.0)

#define TOIDX(x, y) (3*(x + y*WIDTH))

guchar data[3*WIDTH*HEIGHT];
GdkPixbuf *pixbuf;
double window_x0 = -10;
double window_y0 = -10;
double window_x1 = 10;
double window_y1 = 10;
double xscale;
double yscale;
uint16_t nfev = 0;
uint16_t niev = 0;

function function_list[2048];
variable variable_list[128];
expression expression_list[100];
double stack[65536];
double lstack[65536];
char stringbuf[500];
file_data fd;

uint32_t n_expr;
uint32_t n_stack;
expression *top_expr;

clock_t t1, t2;

uint8_t click_state = 0;
double click_x;
double click_y;
int ticker_target, ticker_step;
uint8_t run_ticker;

cairo_matrix_t transform_matrix;


uint32_t eval_func(double t, double *x, double *y, function *func, double *stackpos) {
    variable_list[0].pointer = &t;
    variable_list[0].type = 1<<8;
    uint32_t type = func->oper(func, stackpos);
    if ((type & TYPE_MASK) == TYPE_POINT) {
        *x = stackpos[0];
        *y = stackpos[1];
    } else {
        *y = stackpos[0];
        *x = t;
    }
    nfev++;
    return type;
}

uint32_t eval_inter(double *xi, function *func, double *hstackpos, double *lstackpos) {
    variable_list[0].pointer = xi;
    variable_list[0].type = 1<<8;
    uint32_t type = func->inter(func, hstackpos, lstackpos);
    if ((type & TYPE_MASK) != TYPE_POINT) {
        hstackpos[1] = hstackpos[0];
        hstackpos[0] = xi[1];
        lstackpos[1] = lstackpos[0];
        lstackpos[0] = xi[0];
    }
    niev++;
    return type;
}

void est_radius(double x0, double y0, double x1, double y1, double x2, double y2, double *radius, double *speed) {
    double d01 = (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);
    double d02 = (x0 - x2)*(x0 - x2) + (y0 - y2)*(y0 - y2);
    double d12 = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
    double s = d01 + d02 - d12;
    *radius = sqrt((d01*d02*d12)/(4*d01*d02 - s*s));
    *speed = sqrt(d01);
}

/*void draw_function(function *func, double *stackpos, uint8_t *color) {
    double x, y, xp, yp, xpp, ypp, radius, speed, dt, newdt;
    clock_t t1, t2;
    t1 = clock();
    double t = window_x0;
    nfev = 0;
    eval_func(t, &xpp, &ypp, func, stackpos);
    eval_func(t+DT_DERIV, &xp, &yp, func, stackpos);
    eval_func(t+DT_DERIV*2, &x, &y, func, stackpos);
    est_radius(xpp, ypp, xp, yp, x, y, &radius, &speed);
    dt = fmin(RADIUS_SCALE * radius, MAX_ARC_LENGTH*DT_DERIV / speed);
    eval_func(t+dt, &xp, &yp, func, stackpos);
    
    t = t + dt;
    //draw_line(SCALE_X(xpp), SCALE_Y(ypp), SCALE_X(xp), SCALE_Y(yp));
    LineAA((uint8_t*)data, WIDTH, HEIGHT, SCALE_XK(xpp), SCALE_YK(ypp), SCALE_XK(xp), SCALE_YK(yp), color);
    double minstep_x = (window_x1 - window_x0)/WIDTH;
    double minstep_y = (window_y1 - window_y0)/HEIGHT;
    while (t < window_x1) {
        if (dt < minstep_x) dt = minstep_x;
        eval_func(t+dt, &x, &y, func, stackpos);
        //if (fabs(250*250*((x-xp)*(yp-ypp) - (xp-xpp)*(y-yp))) > 10) {
        //    printf("retrying\n");
        //    dt = dt/2;
        //    continue;
        //}
        t = t + dt;
        //draw_line(SCALE_X(xp), SCALE_Y(yp), SCALE_X(x), SCALE_Y(y));
        LineAA((uint8_t*)data, WIDTH, HEIGHT, SCALE_XK(xp), SCALE_YK(yp), SCALE_XK(x), SCALE_YK(y), color);
        //est_radius(x, y, xp, yp, xpp, ypp, &radius, &speed);
        est_radius(SCALE_XK(x), SCALE_YK(y), SCALE_XK(xp), SCALE_YK(yp), SCALE_XK(xpp), SCALE_YK(ypp), &radius, &speed);
        radius = radius/(1<<16); speed = speed/(1<<16);
        dt = fmin(RADIUS_SCALE * radius, MAX_ARC_LENGTH * dt / speed);
        printf("t: %f -> (%f, %f), radius: %f, arc: %f, dt: %f, det: %f\n", t, x, y, radius, speed, dt, 250*250*((x-xp)*(yp-ypp) - (xp-xpp)*(y-yp)));
        xpp = xp; ypp = yp; xp = x; yp = y;
    }

    t2 = clock();
    g_print("Drawing took %lu, nfev: %d\n", t2-t1, nfev);
}*/

void draw_function_slow(function *func, double *stackpos, uint8_t *color, cairo_t *cr) {
    SET_COLOR(cr, color);
    double x, y, xp, yp, dt, t, t_end;
    if ((eval_func(window_x0, &xp, &yp, func, stackpos) & TYPE_MASK) == TYPE_POINT) {
        t = 0;
        nfev = 0;
        printf("Parametric function detected\n");
        eval_func(0, &xp, &yp, func, stackpos);
        cairo_move_to(cr, SCALE_XK(xp), SCALE_YK(yp));
        eval_func(DT_DERIV, &x, &y, func, stackpos);
        dt = MAX_ARC_LENGTH*DT_DERIV/hypot((x-xp)*xscale, (y-yp)*yscale);
        printf("Initial point (%f, %f), dt is %f\n", xp, yp, dt);
        t += dt;
        while (t < 1) {
            eval_func(t, &x, &y, func, stackpos);
            cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
            dt = MAX_ARC_LENGTH*dt/hypot((x-xp)*xscale, (y-yp)*yscale);
            t += dt;
            xp = x; yp = y;
        }
        printf("Parametric evaluation done in %d steps\n", nfev);
        cairo_stroke(cr);
        return;
    }
    cairo_move_to(cr, SCALE_XK(xp), SCALE_YK(yp));
    double minstep_x = (window_x1 - window_x0)/WIDTH;
    for (t = window_x0+minstep_x; t < window_x1; t+=minstep_x) {
        eval_func(t, &x, &y, func, stackpos);
        if ((x+1 > x) && (y+1 > y) && (xp+1 > xp) && (yp+1 > yp) && (IN_BOUNDS(x, y) || IN_BOUNDS(xp, yp))) cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
        else if ((x+1 > x) && (y+1 > y)) cairo_move_to(cr, SCALE_XK(x), SCALE_YK(y));
        xp = x; yp = y;
    }
    cairo_stroke(cr);
}

void draw_function_constant_ds_rec(function *func, file_data *fd, cairo_t *cr, uint8_t flags, double t_start, double t_end, double xi, double yi, double min_dt) {
    double x, y, xp, yp, dt;
    //printf("t_end %f, t_start %f, min_dt %f\n", t_end, t_start, min_dt);
    if (t_end - t_start < min_dt) {
        //printf("Returning due to minimum dt\n");
        return;
    }
    double *stackpos = fd->stack + fd->n_stack;
    double *lstackpos = fd->lstack;
    //printf("recursive step from %f to %f, flags %02x\n", t_start, t_end, flags);
    if (flags & 0x01) {
        // If the lower point is in-bounds, start from there and go up
        eval_func(t_start, &xp, &yp, func, stackpos);
        cairo_move_to(cr, SCALE_XK(xp), SCALE_YK(yp));
        eval_func(t_start+DT_DERIV, &x, &y, func, stackpos);
        dt = MAX_ARC_LENGTH*DT_DERIV/hypot((x-xp)*xscale, (y-yp)*yscale);
        //printf("Initial point (%f, %f), dt is %f\n", xp, yp, dt);
        t_start += dt;
        while (t_start < t_end) {
            eval_func(t_start, &x, &y, func, stackpos);
            cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
            dt = MAX_ARC_LENGTH*dt/hypot((x-xp)*xscale, (y-yp)*yscale);
            t_start += dt;
            xp = x; yp = y;
            if (!IN_BOUNDS(x, y)) {
                // Leaving the bounds, so iterate on [t, t_end]
                cairo_stroke(cr);
                draw_function_constant_ds_rec(func, fd, cr, (flags & 0x02), t_start, t_end, x, y, min_dt);
                break;
            }
        }
        cairo_stroke(cr);
    } else if (flags & 0x02) {
        // If the upper point is in-bounds, start from there and go down
        eval_func(t_end, &xp, &yp, func, stackpos);
        cairo_move_to(cr, SCALE_XK(xp), SCALE_YK(yp));
        eval_func(t_end-DT_DERIV, &x, &y, func, stackpos);
        dt = MAX_ARC_LENGTH*DT_DERIV/hypot((x-xp)*xscale, (y-yp)*yscale);
        //printf("Initial point (%f, %f), dt is %f\n", xp, yp, dt);
        t_end -= dt;
        while (t_end > t_start) {
            eval_func(t_end, &x, &y, func, stackpos);
            cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
            dt = MAX_ARC_LENGTH*dt/hypot((x-xp)*xscale, (y-yp)*yscale);
            t_end -= dt;
            xp = x; yp = y;
            if (!IN_BOUNDS(x, y)) {
                // Leaving the bounds, so iterate on [t, t_end]
                cairo_stroke(cr);
                draw_function_constant_ds_rec(func, fd, cr, (flags & 0x01), t_start, t_end, x, y, min_dt);
                break;
            }
        }
        cairo_stroke(cr);
    } else if (func->inter) {
        double tempdata[2] = {t_start, t_end};
        eval_inter(tempdata, func, stackpos, lstackpos);
        //printf("Neither point in bounds for [%f, %f] but the interval is ([%f, %f], [%f, %f])\n", t_start, t_end, lstackpos[0], stackpos[0], lstackpos[1], stackpos[1]);
        if ((lstackpos[1] > window_y1) || (stackpos[1] < window_y0) || (lstackpos[0] > window_x1) || (stackpos[0] < window_x0)) return;
        // If neither point is in-bounds, subdivide and iterate
        double t_avg = (t_end + t_start)/2;
        eval_func(t_avg, &x, &y, func, stackpos);
        if (hypot((x - xi)*xscale, (y - yi)*yscale) < MAX_ARC_LENGTH) return;
        //printf("Neither point in bounds on interval %f, %f -> (%f, %f), (%f, %f), %d\n", t_start, t_end, x, y, xi, yi, IN_BOUNDS(x, y));
        if (IN_BOUNDS(x, y)) {
            // [t_start, t_avg]
            draw_function_constant_ds_rec(func, fd, cr, (flags & 0x01) | 0x02, t_start, t_avg, x, y, min_dt);
            // [t_avg, t_end]
            draw_function_constant_ds_rec(func, fd, cr, (flags & 0x02) | 0x01, t_avg, t_end, x, y, min_dt);
        } else {
            draw_function_constant_ds_rec(func, fd, cr, (flags & 0x01), t_start, t_avg, x, y, min_dt);
            draw_function_constant_ds_rec(func, fd, cr, (flags & 0x02), t_avg, t_end, x, y, min_dt);
        }
    } else {
        // If neither point is in-bounds, subdivide and iterate
        double t_avg = (t_end + t_start)/2;
        eval_func(t_avg, &x, &y, func, stackpos);
        if (hypot((x - xi)*xscale, (y - yi)*yscale) < MAX_ARC_LENGTH) return;
        //printf("Neither point in bounds on interval %f, %f -> (%f, %f), (%f, %f), %d\n", t_start, t_end, x, y, xi, yi, IN_BOUNDS(x, y));
        if (IN_BOUNDS(x, y)) {
            // [t_start, t_avg]
            draw_function_constant_ds_rec(func, fd, cr, (flags & 0x01) | 0x02, t_start, t_avg, x, y, min_dt);
            // [t_avg, t_end]
            draw_function_constant_ds_rec(func, fd, cr, (flags & 0x02) | 0x01, t_avg, t_end, x, y, min_dt);
        } else {
            draw_function_constant_ds_rec(func, fd, cr, (flags & 0x01), t_start, t_avg, x, y, min_dt);
            draw_function_constant_ds_rec(func, fd, cr, (flags & 0x02), t_avg, t_end, x, y, min_dt);
        }
    }
}

void draw_function_constant_ds(function *func, file_data *fd, uint8_t *color, cairo_t *cr) {
    SET_COLOR(cr, color);
    double x, y, xp, yp, dt, t, t_end;
    double *stackpos = fd->stack + fd->n_stack;
    t = 0;
    t_end = 1;
    double min_dt = MIN_DT;
    if ((eval_func(0, &xp, &yp, func, stackpos) & TYPE_MASK) != TYPE_POINT) {
        t = window_x0;
        t_end = window_x1;
        min_dt = MAX_ARC_LENGTH/xscale;
    }
    nfev = 0;
    niev = 0;
    uint8_t flags = 0;
    eval_func(t, &xp, &yp, func, stackpos);
    if (IN_BOUNDS(xp, yp)) flags |= 0x01;
    eval_func(t_end, &x, &y, func, stackpos);
    if (IN_BOUNDS(x, y)) flags |= 0x02;
    draw_function_constant_ds_rec(func, fd, cr, flags, t, t_end, xp, yp, min_dt);
}

void draw_implicit_rec(function *func, file_data *fd, cairo_t *cr, double *area, int divisions) {
    fd->variable_list[0].pointer = area;
    fd->variable_list[0].type = 1<<8;
    fd->variable_list[1].pointer = area+2;
    fd->variable_list[1].type = 1<<8;
    double *stackpos = fd->stack + fd->n_stack;
    double *lstackpos = fd->lstack;
    func->inter(func, stackpos, lstackpos);
    niev++;
    if ((lstackpos[0] > 0) || (stackpos[0] < 0)) {
        // No contours can exist in the given area
        //printf("    No contours found in ([%f, %f], [%f, %f])\n", area[0], area[1], area[2], area[3]);
        //cairo_rectangle(cr, SCALE_XK(area[0]), SCALE_YK(area[3]), xscale*(area[1] - area[0]), yscale*(area[3] - area[2]));
        //cairo_stroke(cr);
        return;
    }
    if (divisions == IMPLICIT_MAXDEPTH) {
        // Maximum depth reached
        //printf("    maximum depth reached on ([%f, %f], [%f, %f]) --> [%f, %f]\n", area[0], area[1], area[2], area[3], lstackpos[0], stackpos[0]);
        function *arg1 = func->first_arg;
        function *arg2 = arg1->next_arg;
        double x0 = area[0], x1 = area[1], y0 = area[2], y1 = area[3];
        fd->variable_list[0].pointer = &x0;
        fd->variable_list[1].pointer = &y0;
        arg1->oper(arg1, stackpos);
        double e00 = stackpos[0];
        arg2->oper(arg2, stackpos);
        e00 -= stackpos[0];
        fd->variable_list[0].pointer = &x1;
        fd->variable_list[1].pointer = &y0;
        arg1->oper(arg1, stackpos);
        double e10 = stackpos[0];
        arg2->oper(arg2, stackpos);
        e10 -= stackpos[0];
        fd->variable_list[0].pointer = &x1;
        fd->variable_list[1].pointer = &y1;
        arg1->oper(arg1, stackpos);
        double e11 = stackpos[0];
        arg2->oper(arg2, stackpos);
        e11 -= stackpos[0];
        fd->variable_list[0].pointer = &x0;
        fd->variable_list[1].pointer = &y1;
        arg1->oper(arg1, stackpos);
        double e01 = stackpos[0];
        arg2->oper(arg2, stackpos);
        e01 -= stackpos[0];
        nfev += 4;
        double npos, epos, spos, wpos;
        uint8_t edges = 0;
        if ((e00 != e10) && (((e00 <= 0) && (e10 >= 0)) || ((e00 >= 0) && (e10 <= 0)))) {
            // south edge
            edges |= 0x4;
            spos = x0 - e00*(x1 - x0)/(e10 - e00);
        }
        if ((e10 != e11) && (((e10 <= 0) && (e11 >= 0)) || ((e10 >= 0) && (e11 <= 0)))) {
            // east edge
            edges |= 0x2;
            epos = y0 - e10*(y1 - y0)/(e11 - e10);
        }
        if ((e11 != e01) && (((e11 <= 0) && (e01 >= 0)) || ((e11 >= 0) && (e01 <= 0)))) {
            // north edge
            edges |= 0x1;
            npos = x0 - e01*(x0 - x1)/(e01 - e11);
        }
        if ((e01 != e00) && (((e01 <= 0) && (e00 >= 0)) || ((e01 >= 0) && (e00 <= 0)))) {
            // west edge
            edges |= 0x8;
            wpos = y0 - e00*(y0 - y1)/(e00 - e01);
        }
        // If two edges are selected, then draw the line between those two edges. If three 
        // edges are selected, then one of the corners must be zero so the contour must go
        // through one of the corners. If four edges are selected, there are two contour
        // lines, each going through opposite sides of the area.
        if ((edges & 0x5) == 0x5) {
            cairo_move_to(cr, SCALE_XK(spos), SCALE_YK(y0));
            cairo_line_to(cr, SCALE_XK(npos), SCALE_YK(y1));
            cairo_stroke(cr);
            edges &= 0xa;
        }
        if ((edges & 0xa) == 0xa) {
            cairo_move_to(cr, SCALE_XK(x0), SCALE_YK(wpos));
            cairo_line_to(cr, SCALE_XK(x1), SCALE_YK(epos));
            cairo_stroke(cr);
            edges &= 0x5;
        }
        if ((edges & 0x3) == 0x3) {
            cairo_move_to(cr, SCALE_XK(x1), SCALE_YK(epos));
            cairo_line_to(cr, SCALE_XK(npos), SCALE_YK(y1));
            cairo_stroke(cr);
        } else if ((edges & 0x6) == 0x6) {
            cairo_move_to(cr, SCALE_XK(x1), SCALE_YK(epos));
            cairo_line_to(cr, SCALE_XK(spos), SCALE_YK(y0));
            cairo_stroke(cr);
        } else if ((edges & 0xc) == 0xc) {
            cairo_move_to(cr, SCALE_XK(x0), SCALE_YK(wpos));
            cairo_line_to(cr, SCALE_XK(spos), SCALE_YK(y0));
            cairo_stroke(cr);
        } else if ((edges & 0x9) == 0x9) {
            cairo_move_to(cr, SCALE_XK(x0), SCALE_YK(wpos));
            cairo_line_to(cr, SCALE_XK(npos), SCALE_YK(y1));
            cairo_stroke(cr);
        }

    } else {
        // Subdivide
        double x0 = area[0], x1 = area[1], y0 = area[2], y1 = area[3];
        double xm = (x0 + x1)/2, ym = (y0 + y1)/2;
        //printf("subdividing ([%f, %f], [%f, %f]), %f, %f\n", x0, x1, y0, y1, xm, ym);
        area[1] = xm; area[3] = ym;
        draw_implicit_rec(func, fd, cr, area, divisions+1);
        area[0] = xm; area[1] = x1; area[2] = y0; area[3] = ym;
        draw_implicit_rec(func, fd, cr, area, divisions+1);
        area[0] = x0; area[1] = xm; area[2] = ym; area[3] = y1;
        draw_implicit_rec(func, fd, cr, area, divisions+1);
        area[0] = xm; area[1] = x1; area[2] = ym; area[3] = y1;
        draw_implicit_rec(func, fd, cr, area, divisions+1);
    }
}

void draw_implicit(function *func, file_data *fd, uint8_t *color, cairo_t *cr) {
    SET_COLOR(cr, color);
    double temp[4] = {window_x0, window_x1, window_y0, window_y1};
    draw_implicit_rec(func, fd, cr, temp, 0);
}

static gboolean button_press_callback (GtkWidget *event_box, GdkEventButton *event, gpointer data) {
    g_print ("Clicked at %f, %f\n", event->x, event->y);
    click_x = event->x;
    click_y = event->y;
    click_state = 1;
    return TRUE;
}

static gboolean button_release_callback (GtkWidget *event_box, GdkEventButton *event, gpointer data) {
    click_x = event->x;
    click_y = event->y;
    click_state = 0;
    return TRUE;
}

gboolean redraw_all(GtkWidget *widget, cairo_t *cr, gpointer data_pointer) {
    file_data *fd = (file_data*)(data_pointer);
    cairo_set_source_rgb(cr, COLOR_MAJOR_GRID/256.0, COLOR_MAJOR_GRID/256.0, COLOR_MAJOR_GRID/256.0);
    cairo_set_line_width(cr, 1.0);
    
    xscale = 1.0*WIDTH/(window_x1 - window_x0);
    yscale = 1.0*HEIGHT/(window_y1 - window_y0);
    t1 = clock();
    
    int16_t logx = round(3*log10(TICK_SIZE / xscale));
    int8_t basex = logx%3;
    basex = (basex<0 ? 3+basex : basex);
    int16_t decadex = (logx - basex)/3;
    basex++;
    if (basex == 3) basex = 5;
    clock_t t3, t4;
    double ticksize = basex;
    int64_t xk;
    while (decadex > 0) {decadex--; ticksize *= 10;}
    while (decadex < 0) {decadex++; ticksize /= 10;}
    double x0_scaled = ticksize*((int64_t)(window_x0/ticksize));
    double x1_scaled = ticksize*((int64_t)(window_x1/ticksize));
    for (double tick=x0_scaled; tick <= x1_scaled; tick+=ticksize) {
        cairo_move_to(cr, SCALE_XK(tick), 0);
        cairo_line_to(cr, SCALE_XK(tick), HEIGHT);
    }
    
    int16_t logy = round(3*log10(TICK_SIZE / yscale));
    int8_t basey = logy%3;
    basey = (basey<0 ? 3+basey : basey);
    int16_t decadey = (logy - basey)/3;
    basey++;
    if (basey == 3) basey = 5;
    ticksize = basey;
    while (decadey > 0) {decadey--; ticksize *= 10;}
    while (decadey < 0) {decadey++; ticksize /= 10;}
    double y0_scaled = ticksize*((int64_t)(window_y0/ticksize));
    double y1_scaled = ticksize*((int64_t)(window_y1/ticksize));
    for (double tick=y0_scaled; tick <= y1_scaled; tick+=ticksize) {
        cairo_move_to(cr, 0, SCALE_YK(tick));
        cairo_line_to(cr, WIDTH, SCALE_YK(tick));
    }

    cairo_stroke(cr);
    // Draw axes gridlines
    cairo_set_source_rgb(cr, COLOR_AXES/256.0, COLOR_AXES/256.0, COLOR_AXES/256.0);
    if ((window_x0 <= 0) && (0 < window_x1)) {
        cairo_move_to(cr, SCALE_XK(0), 0);
        cairo_line_to(cr, SCALE_XK(0), HEIGHT);
    }
    if ((window_y0 < 0) && (0 <= window_y1)) {
        cairo_move_to(cr, 0, SCALE_YK(0));
        cairo_line_to(cr, WIDTH, SCALE_YK(0));
    }
    cairo_stroke(cr);


    double pt_x, pt_y, pt_x1, pt_y1;
    expression *expression_list = fd->expression_list, *expr;
    double *stack = fd->stack;
    uint32_t n_stack = fd->n_stack;
    for (int i=0; i < fd->n_expr; i++) {
        expr = expression_list+i;
        if (expr->flags & EXPRESSION_PLOTTABLE) {
            if ((expr->flags & EXPRESSION_FIXED) && ((expr->value_type & TYPE_MASK) == TYPE_POINT)) {
#ifdef DEBUG_PLOT
                t3 = clock();
#endif
                int len = (expr->value_type) >> 8;
                double *ptr = expr->value;
                SET_COLOR(cr, expr->color);
                for (int p=0; p < len; p+=2) {
                    pt_x = SCALE_XK(ptr[p]);
                    pt_y = SCALE_YK(ptr[p+1]);
                    cairo_arc(cr, pt_x, pt_y, POINT_SIZE, 0, 2*G_PI);
                    cairo_fill(cr);
                }
#ifdef DEBUG_PLOT
                t4 = clock();
                printf("Plotted point expression %p (%d) in %luus\n", expression_list+i, i+1, t4-t3);
#endif
            } else if ((expr->flags & EXPRESSION_FIXED) && ((expr->value_type & TYPE_MASK) == TYPE_POLYGON)) {
#ifdef DEBUG_PLOT
                t3 = clock();
#endif
                int len = (expr->value_type) >> 8;
                double *ptr = expr->value;
                double x_temp, y_temp;
                if (i != 19) SET_COLOR(cr, expr->color);
                uint32_t color_pos = 0;
                for (int p=0; p < len; p += 2*MAX_POLYGON_SIZE) {
                    pt_x = SCALE_XK(ptr[p]);
                    pt_y = SCALE_YK(ptr[p+1]);
                    //if (i == 19) SET_COLOR(cr, (expression_list[20].value+color_pos));
                    color_pos += 3;
                    cairo_move_to(cr, pt_x, pt_y);
                    for (int k=1; k < MAX_POLYGON_SIZE; k++) {
                        x_temp = ptr[p+2*k]; y_temp = ptr[p+2*k+1];
                        if ((x_temp+1 > x_temp) && (y_temp+1 > y_temp)) {
                            pt_x1 = SCALE_XK(x_temp);
                            pt_y1 = SCALE_YK(y_temp);
                            cairo_line_to(cr, pt_x1, pt_y1);
                            pt_x = pt_x1;
                            pt_y = pt_y1;
                        }
                    }
                    cairo_close_path(cr);
                    cairo_stroke_preserve(cr);
                    cairo_fill(cr);
                }
#ifdef DEBUG_PLOT
                t4 = clock();
                printf("Plotted polygon expression %p (%d) in %luus\n", expression_list+i, i+1, t4-t3);
#endif
            } else if (expr->func->oper == func_equals) {
                if (!(expr->func->inter)) {
                    printf("ERROR: interval function needed for implicit plotting\n");
                    exit(EXIT_FAILURE);
                }
#ifdef DEBUG_PLOT
                t3 = clock();
                nfev = 0;
                niev = 0;
#endif
                draw_implicit(expression_list[i].func, fd, expr->color, cr);
#ifdef DEBUG_PLOT
                t4 = clock();
                printf("Plotted implicit expression %p (%d) in %luus, nfev: %d, niev: %d\n", expression_list+i, i+1, t4-t3, nfev, niev);
#endif
            } else {
#ifdef DEBUG_PLOT
                t3 = clock();
                nfev = 0;
                niev = 0;
#endif
                draw_function_constant_ds(expression_list[i].func, fd, expr->color, cr);
#ifdef DEBUG_PLOT
                t4 = clock();
                printf("Plotted expression %p (%d) in %luus, nfev: %d, niev: %d, interval function %p\n", expression_list+i, i+1, t4-t3, nfev, niev, expression_list[i].func->inter);
#endif
            }
        }
    }
    t2 = clock();
    g_print("Redraw took %luus, %d expressions, bounds %f %f %f %f\n", t2-t1, n_expr, window_x0, window_y0, window_x1, window_y1);
    return FALSE;
}

static gboolean scroll_callback (GtkWidget *event_box, GdkEventScroll *event, gpointer data_pointer) {
    //g_print ("Scrolled at %f, %f, direction %d, data_pointer: %p\n", event->x, event->y, event->direction, data_pointer);
    double center_x = ((int32_t)(event->x))/xscale + window_x0;
    double center_y = window_y1 - ((int32_t)(event->y))/yscale;
    //g_print("center %f, %f, bounds %f, %f, %f, %f\n", center_x, center_x, window_x0, window_y0, window_x1, window_y1);
    uint8_t redraw = 0;
    // 1 = scroll out, 0 = scroll in
    if (event->direction == 0) {
        window_x0 = center_x - (center_x - window_x0)/SCROLL_AMOUNT;
        window_x1 = center_x - (center_x - window_x1)/SCROLL_AMOUNT;
        window_y0 = center_y - (center_y - window_y0)/SCROLL_AMOUNT;
        window_y1 = center_y - (center_y - window_y1)/SCROLL_AMOUNT;
        redraw = 1;
    } else if (event->direction == 1) {
        window_x0 = center_x - (center_x - window_x0)*SCROLL_AMOUNT;
        window_x1 = center_x - (center_x - window_x1)*SCROLL_AMOUNT;
        window_y0 = center_y - (center_y - window_y0)*SCROLL_AMOUNT;
        window_y1 = center_y - (center_y - window_y1)*SCROLL_AMOUNT;
        redraw = 1;
    }

    if (redraw) {
        gtk_widget_queue_draw(((file_data*)(data_pointer))->drawing_area);
    }
    return TRUE;
}

static gboolean motion_callback(GtkWidget *event_box, GdkEventMotion *event, gpointer data_pointer) {
    if (click_state) {
        double dx = (event->x - click_x)/xscale;
        double dy = (click_y - event->y)/yscale;
        window_x0 -= dx; window_x1 -= dx; window_y0 -= dy; window_y1 -= dy;
        click_x = event->x;
        click_y = event->y;
        uint8_t color[4] = {255, 255, 255, 0};
        gtk_widget_queue_draw(((file_data*)(data_pointer))->drawing_area);
    }
    return TRUE;
}

static gboolean timeout_callback(gpointer data_pointer) {
    printf("timeout_callback\n");
    clock_t t3 = clock();
    file_data *fd = (file_data*)(data_pointer);
    fd->expression_list[ticker_target].func->oper(fd->expression_list[ticker_target].func, (fd->stack)+(fd->n_stack));
    expression *expr = top_expr;
    expression *from = NULL;
    while (expr) {
        if ((expr->var) && (expr->var->new_pointer)) {
            printf("expression %p (offset %ld) has changed, expr->var %p\n", expr, expr - (fd->expression_list) + 1, expr->var);
            if (!from) from = expr;
            expr->var->pointer = expr->var->new_pointer;
            expr->var->type = expr->var->new_type;
            expr->var->new_pointer = NULL;
            // If we assign to a variable, we must unlink the function block,
            // or else the value will be overwritten
            expr->func = NULL;
        }
        expr = expr->next_expr;
    }
    if (from) {
        printf("evaluating from %p\n", from);
        evaluate_from(fd, from);
        clock_t t4 = clock();
        printf("Total evaluation time: %luus\n", t4-t3);
        gtk_widget_queue_draw(fd->drawing_area);
    }
    if (run_ticker) return TRUE;
    return FALSE;
}

static gboolean keypress_callback(GtkWidget *widget, GdkEventKey *event, gpointer data_pointer) {
    printf("Event is %d\n", event->keyval);
    if (event->keyval == 's') {
        run_ticker = 1;
        if (ticker_target >= 0) g_timeout_add(ticker_step, timeout_callback, data_pointer);
    } else if (event->keyval == 'e') {
        run_ticker = 0;
    } else if (event->keyval == 'i') {
        treeview_activate(&fd);
        printf("inspecting\n");
    }
    return TRUE;
}

/*void eval_func(double t, double *x, double *y) {
    double e = exp(t);
    *x = cos(e)/e;
    *y = sin(e)/e;
}*/

/*void eval_func(double t, double *x, double *y) {
    *x = 2*cos(t) + 5*cos(2*t/3);
    *y = 2*sin(t) + 5*sin(2*t/3);
}*/

/*void eval_func(double t, double *x, double *y) {
    nfev++;
    *x = 2*(sin(5*t) - sin(3*t));
    *y = 2*sin(cos(7*t) + cos(t));
}*/


static void activate (GtkApplication *app, gpointer user_data) {
    GtkWidget *window;
    GtkWidget *event_box;
    GtkWidget *image;
    file_data *fd = (file_data*)(user_data);

    window = gtk_application_window_new (app);
    event_box = gtk_event_box_new();
    gtk_container_add(GTK_CONTAINER(window), event_box);

    gtk_widget_add_events(event_box, GDK_BUTTON_PRESS_MASK);
    gtk_widget_add_events(event_box, GDK_SCROLL_MASK);
    gtk_widget_add_events(event_box, GDK_POINTER_MOTION_MASK);
    gtk_widget_add_events(event_box, GDK_BUTTON_RELEASE_MASK);

    GtkWidget *drawing_area = gtk_drawing_area_new();
    gtk_widget_set_size_request(drawing_area, WIDTH, HEIGHT);
    fd->drawing_area = drawing_area;
    g_signal_connect(G_OBJECT(event_box), "button_press_event", G_CALLBACK(button_press_callback), fd);
    g_signal_connect(G_OBJECT(event_box), "button_release_event", G_CALLBACK(button_release_callback), fd);
    g_signal_connect(G_OBJECT(event_box), "scroll_event", G_CALLBACK(scroll_callback), fd);
    g_signal_connect(G_OBJECT(event_box), "motion-notify-event", G_CALLBACK(motion_callback), fd);
    g_signal_connect(G_OBJECT(drawing_area), "draw", G_CALLBACK(redraw_all), fd);
    g_signal_connect(window, "key-press-event", G_CALLBACK(keypress_callback), fd);
    gtk_container_add(GTK_CONTAINER(event_box), drawing_area);

    gtk_widget_show_all(window);
}

int main (int argc, char **argv) {
    memset(variable_list, 0, 10*sizeof(variable));
    memset(expression_list, 0, 10*sizeof(expression));
    printf("First stack positions %p, %p\n", stack, stack+1);
    variable_list[0] = new_variable("x", 0, VARIABLE_IN_SCOPE, NULL);
    variable_list[1] = new_variable("y", 0, VARIABLE_IN_SCOPE, NULL);
    double pi = M_PI;
    variable_list[2] = new_variable("\\pi", 1<<8, VARIABLE_IN_SCOPE, &pi);

    fd.expression_list = expression_list;
    fd.variable_list = variable_list;
    fd.function_list = function_list;
    fd.stack = stack;
    fd.lstack = lstack;
    fd.n_expr = 0;
    fd.n_var = 3;
    fd.n_func = 0;
    fd.n_stack = 0;
    
    uint32_t n_func = 0;
    uint32_t n_var = 0;
    load_file(argv[1], &fd);
    top_expr = parse_file(&fd, stringbuf);
    printf("Parsing completed. %d function blocks, %d variables, %d expressions\n", fd.n_func, fd.n_var, fd.n_expr);
    ticker_target = (int)parse_double(argv[2]);
    ticker_step = (int)parse_double(argv[3]);
    printf("Expression %d will be evaluated every %d milliseconds\n", ticker_target, ticker_step);

    uint32_t type;

    GtkApplication *app;
    int status;


    app = gtk_application_new("org.gtk.example", G_APPLICATION_DEFAULT_FLAGS);
    g_signal_connect(app, "activate", G_CALLBACK (activate), &fd);
    for (int i=4; i < argc; i++) argv[i-3] = argv[i];
    status = g_application_run(G_APPLICATION (app), argc-3, argv);
    g_object_unref (app);

    return status;
}
