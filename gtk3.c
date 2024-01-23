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
#define IN_BOUNDS(x, y) ((window_x0 <= x) && (x <= window_x1) && (window_y0 <= y) && (y <= window_y1) && (x == x) && (y == y))
#define SET_COLOR(cr, color) cairo_set_source_rgb(cr, color[0]/256.0, color[1]/256.0, color[2]/256.0)
#define SET_COLOR_OPACITY(cr, color, opacity) cairo_set_source_rgba(cr, color[0]/256.0, color[1]/256.0, color[2]/256.0, opacity)

#define TOIDX(x, y) (3*(x + y*WIDTH))

#define GET_STEP(type) (step_table[(type) & TYPE_MASK])
const static uint32_t step_table[8] = {1, 2, 3, 1, 1, MAX_POLYGON_SIZE*2, 0, 0};

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
variable variable_list[512];
expression expression_list[100];
expression *deptable[300];
double stack[65536];
double lstack[65536];
char stringbuf[600];
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

uint32_t eval_index = 0;

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

uint32_t eval_inter_dep(expression *expr, double *hstackpos, double *lstackpos) {
    //if (expr->var) printf("Evaluating interval for %p (%s), non-fixed %d, num %d\n", expr, expr->var->name, expr->num_nonfixed_dependencies, expr->num_dependencies);
    uint32_t type;
    uint32_t len;
    expression *dep;
    if (expr->num_nonfixed_dependencies != 0) {
        for (int i=0; i < expr->num_dependencies; i++) {
            dep = expr->dependencies[i];
            if (dep->flags & EXPRESSION_FIXED) continue;
            type = eval_inter_dep(dep, hstackpos, lstackpos);
            if (dep->var->flags & VARIABLE_FUNCTION) continue;
            len = type>>8;
            if (!(dep->var->pointer) || (dep->var->type != type)) dep->var->pointer = realloc(dep->var->pointer, 2*len*sizeof(double));
            memcpy(dep->var->pointer, lstackpos, len*sizeof(double));
            memcpy(dep->var->pointer+len, hstackpos, len*sizeof(double));
            //printf("interval is %08x\n", type);
            //print_object(type, lstackpos); printf(" to "); print_object(type, hstackpos); printf("\n");
            dep->var->type = type;
            dep->var->flags |= VARIABLE_INTERVAL;
        }
    }
    if ((expr->var) && (expr->var->flags & VARIABLE_FUNCTION)) return 0;
    //printf("calling interval function %p from block %p\n", expr->func->inter, expr->func);
    type = expr->func->inter(expr->func, hstackpos, lstackpos);
    //if (expr->num_nonfixed_dependencies) printf("result is %08x: [%f, %f]\n", type, lstackpos[0], hstackpos[0]);
    return type;
}

uint32_t eval_inter_2d(double *x, double *y, expression *expr, double *hstackpos, double *lstackpos) {
    variable_list[0].pointer = x;
    variable_list[0].type = 1<<8;
    variable_list[1].pointer = y;
    variable_list[1].type = 1<<8;
    niev++;
    return eval_inter_dep(expr, hstackpos, lstackpos);
}

uint32_t eval_func_dep(expression *expr, double *stackpos) {
    uint32_t type;
    uint32_t len;
    expression *dep;
    if (expr->num_nonfixed_dependencies != 0) {
        for (int i=0; i < expr->num_dependencies; i++) {
            dep = expr->dependencies[i];
            if (dep->flags & EXPRESSION_FIXED) continue;
            type = eval_func_dep(dep, stackpos);
            if (dep->var->flags & VARIABLE_FUNCTION) continue;
            len = type>>8;
            if (!(dep->var->pointer) || (dep->var->type != type)) dep->var->pointer = realloc(dep->var->pointer, len*sizeof(double));
            memcpy(dep->var->pointer, stackpos, len*sizeof(double));
            dep->var->type = type;
        }
    }
    if ((expr->var) && (expr->var->flags & VARIABLE_FUNCTION)) return 0;
    //printf("calling interval function %p from block %p\n", expr->func->inter, expr->func);
    type = expr->func->oper(expr->func, stackpos);
    //if (expr->num_nonfixed_dependencies) printf("result is %08x: [%f, %f]\n", type, lstackpos[0], hstackpos[0]);
    return type;
}

uint32_t eval_func_2d(double *x, double *y, expression *expr, double *stackpos) {
    variable_list[0].pointer = x;
    variable_list[0].type = 1<<8;
    variable_list[1].pointer = y;
    variable_list[1].type = 1<<8;
    nfev++;
    return eval_func_dep(expr, stackpos);
}

uint32_t eval_func(double t, double *x, double *y, expression *expr, double *stackpos) {
    variable_list[0].pointer = &t;
    variable_list[0].type = 1<<8;
    //uint32_t type = expr->func->oper(expr->func, stackpos);
    uint32_t type = eval_func_dep(expr, stackpos);
    if ((type & TYPE_MASK) == TYPE_POINT) {
        *x = stackpos[2*eval_index];
        *y = stackpos[2*eval_index+1];
    } else {
        *y = stackpos[eval_index];
        *x = t;
    }
    if (expr->cache_size != -1) {
        if (expr->cache_size % CACHE_BLOCK_SIZE == 0) {
            expr->value = realloc(expr->value, (expr->cache_size + CACHE_BLOCK_SIZE) * sizeof(double));
        }
        expr->value[expr->cache_size] = *x;
        expr->value[expr->cache_size+1] = *y;
        expr->cache_size += 2;
    }
    nfev++;
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

void find_extremum(expression *expr, double *stackpos, double x, double y, double xp, double yp, double xpp, double ypp) {
    double u1, v1, u2, v2, ex, temp, val, m1, m2;
    // We assume that x > xp > xpp
    for (int i=0; i < 10; i++) {
        u1 = x - xp; v1 = y - yp;
        u2 = xpp - xp; v2 = ypp - yp;
        temp = fabs(yp) * PLOT_EXTREMUM_EPS;
        if ((fabs(v1) < temp) && (fabs(v2) < temp)) break;
        ex = (u2*u2*v1 - u1*u1*v2)/(2*(u2*v1 - u1*v2));
        if (ex > 0) {
            xpp = xp;
            ypp = yp;
            xp += ex;
            eval_func(xp, &temp, &yp, expr, stackpos);
        } else if (ex < 0) {
            x = xp;
            y = yp;
            xp += ex;
            eval_func(xp, &temp, &yp, expr, stackpos);
        } else {
            x = (x + xp)/2;
            eval_func(x, &temp, &y, expr, stackpos);
            xpp = (xpp + xp)/2;
            eval_func(xpp, &temp, &ypp, expr, stackpos);
        }
        printf("    iterated extremum at %f, value %f\n", xp, yp);
    }
}

void find_discontinuity(expression *expr, double *stackpos, double t, double x, double y, double tp, double xp, double yp, cairo_t *cr) {
    // Assume that the function is monotonic on the interval [t, tp]
    // If we subdivide the interval
    double tc, xc, yc, d, dp;
    double origdist = PLOT_CONTINUITY_THRESHOLD*hypot((x-xp)*xscale, (y-yp)*yscale);
    for (int i=0; i < PLOT_DISCONTINUITY_MAXITER; i++) {
        if (hypot((x-xp)*xscale, (y-yp)*yscale) <= origdist) {
            cairo_line_to(cr, SCALE_XK(xp), SCALE_YK(yp));
            cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
            return;
        }
        tc = (t + tp)/2;
        eval_func(tc, &xc, &yc, expr, stackpos);
        d = hypot((xc-x)*xscale, (yc-y)*yscale);
        dp = hypot((xc-xp)*xscale, (yc-yp)*yscale);
        if (d > dp) {
            // tc is close to tp so the discontinuity is probably in
            // the range [t, tc]
            tp = tc; xp = xc; yp = yc;
        } else {
            // tc is close to t so the discontinuity is probably in
            // the range [tc, tp]
            t = tc; x = xc; y = yc;
        }
        //printf("    iterated discontinuity %f -> (%f, %f) and %f -> (%f, %f)\n", tp, xp, yp, t, x, y);
    }
    cairo_line_to(cr, SCALE_XK(xp), SCALE_YK(yp));
    cairo_move_to(cr, SCALE_XK(x), SCALE_YK(y));
}

void draw_function_constant_ds_rec(expression *expr, file_data *fd, cairo_t *cr, uint8_t flags, double t_start, double t_end, double xi, double yi, double min_dt) {
    double x, y, xp, yp, xpp, ypp, tp, tpp, dt, ds, dsp;
    printf("t_end %f, t_start %f, min_dt %f, flags %02x, interval %p\n", t_end, t_start, min_dt, flags, expr->func->inter);
    if (t_end - t_start < min_dt) {
        //printf("Returning due to minimum dt\n");
        return;
    }
    double *stackpos = fd->stack + fd->n_stack;
    double *lstackpos = fd->lstack;
    //printf("recursive step from %f to %f, flags %02x\n", t_start, t_end, flags);
    if (flags & 0x01) {
        // If the lower point is in-bounds, start from there and go up
        eval_func(t_start, &xp, &yp, expr, stackpos);
        cairo_move_to(cr, SCALE_XK(xp), SCALE_YK(yp));
        eval_func(t_start+PLOT_DT_DERIV, &x, &y, expr, stackpos);
        ds = hypot((x-xp)*xscale, (y-yp)*yscale);
        if (ds > PLOT_MAX_ARC_LENGTH*PLOT_DISCONTINUITY_THRESHOLD) {
            //printf("Possible discontinuity between %f and %f: (%f, %f) and (%f, %f)\n", t_start, tp, x, y, xp, yp);
            find_discontinuity(expr, stackpos, t_start, x, y, tp, xp, yp, cr);
        }
        dt = PLOT_MAX_ARC_LENGTH*PLOT_DT_DERIV/ds;
        //printf("Initial point (%f, %f), dt is %e, dx is %e\n", xp, yp, dt, (x-xp));
        tp = t_start;
        t_start += dt;
        uint32_t step = 0;
        while (t_start < t_end) {
            eval_func(t_start, &x, &y, expr, stackpos);
            if ((x == x) && (y == y)) {
                ds = hypot((x-xp)*xscale, (y-yp)*yscale);
                dt = PLOT_MAX_ARC_LENGTH*dt/ds;
                /*if ((yp > ypp) && (yp > y) && (step >= 2)) {
                    printf("Possible maximum at %f -> (%f, %f)\n", tp, xp, yp);
                    find_extremum(func, stackpos, t_start, y, tp, xp, tpp, ypp);
                } else if ((yp < ypp) && (yp < y) && (step >= 2)) {
                    printf("Possible minimum at %f -> (%f, %f)\n", tp, xp, yp);
                }*/
                if (ds > PLOT_MAX_ARC_LENGTH*PLOT_DISCONTINUITY_THRESHOLD) {
                    //printf("Possible discontinuity between %f and %f: (%f, %f) and (%f, %f)\n", t_start, tp, x, y, xp, yp);
                    find_discontinuity(expr, stackpos, t_start, x, y, tp, xp, yp, cr);
                }
                cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
                step++;
                xpp = xp; ypp = yp; xp = x; yp = y; tpp = tp; tp = t_start; dsp = ds;
                t_start += dt;
            }
            if (!IN_BOUNDS(x, y)) {
                // Leaving the bounds, so iterate on [t, t_end]
                cairo_stroke(cr);
                draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x02), t_start, t_end, x, y, min_dt);
                break;
            }
        }
        cairo_stroke(cr);
    } else if (flags & 0x02) {
        // If the upper point is in-bounds, start from there and go down
        eval_func(t_end, &xp, &yp, expr, stackpos);
        cairo_move_to(cr, SCALE_XK(xp), SCALE_YK(yp));
        eval_func(t_end-PLOT_DT_DERIV, &x, &y, expr, stackpos);
        ds = hypot((x-xp)*xscale, (y-yp)*yscale);
        dt = PLOT_MAX_ARC_LENGTH*PLOT_DT_DERIV/ds;
        if (ds > PLOT_MAX_ARC_LENGTH*PLOT_DISCONTINUITY_THRESHOLD) {
            //printf("Possible discontinuity between %f and %f: (%f, %f) and (%f, %f)\n", t_start, tp, x, y, xp, yp);
            find_discontinuity(expr, stackpos, t_start, x, y, tp, xp, yp, cr);
        }
        //printf("Initial point (%f, %f), dt is %f\n", xp, yp, dt);
        t_end -= dt;
        while (t_end > t_start) {
            eval_func(t_end, &x, &y, expr, stackpos);
            if ((x == x) && (y == y)) {
                ds = hypot((x-xp)*xscale, (y-yp)*yscale);
                dt = PLOT_MAX_ARC_LENGTH*dt/ds;
                if (ds > PLOT_MAX_ARC_LENGTH*PLOT_DISCONTINUITY_THRESHOLD) {
                    //printf("Possible discontinuity between %f and %f: (%f, %f) and (%f, %f)\n", t_end, tp, x, y, xp, yp);
                    find_discontinuity(expr, stackpos, t_end, x, y, tp, xp, yp, cr);
                }
                cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
                xp = x; yp = y; tp = t_end; dsp = ds;
                t_end -= dt;
            }
            if (!IN_BOUNDS(x, y)) {
                // Leaving the bounds, so iterate on [t, t_end]
                cairo_stroke(cr);
                draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x01), t_start, t_end, x, y, min_dt);
                break;
            }
        }
        cairo_stroke(cr);
    } else if (expr->func->inter) {
        double tempdata[2] = {t_start, t_end};
        eval_inter(tempdata, expr->func, stackpos, lstackpos);
        //printf("Neither point in bounds for [%f, %f] but the interval is ([%f, %f], [%f, %f])\n", t_start, t_end, lstackpos[0], stackpos[0], lstackpos[1], stackpos[1]);
        if ((lstackpos[1] > window_y1) || (stackpos[1] < window_y0) || (lstackpos[0] > window_x1) || (stackpos[0] < window_x0) || ((stackpos[0] != stackpos[0]) && (lstackpos[0] != lstackpos[0])) || ((stackpos[1] != stackpos[1]) && (lstackpos[1] != lstackpos[1]))) return;
        // If neither point is in-bounds, subdivide and iterate
        double t_avg = (t_end + t_start)/2;
        eval_func(t_avg, &x, &y, expr, stackpos);
        if (hypot((x - xi)*xscale, (y - yi)*yscale) < PLOT_MAX_ARC_LENGTH) return;
        //printf("Neither point in bounds on interval %f, %f -> (%f, %f), (%f, %f), %d\n", t_start, t_end, x, y, xi, yi, IN_BOUNDS(x, y));
        if (IN_BOUNDS(x, y)) {
            // [t_start, t_avg]
            draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x01) | 0x02, t_start, t_avg, x, y, min_dt);
            // [t_avg, t_end]
            draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x02) | 0x01, t_avg, t_end, x, y, min_dt);
        } else {
            draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x01), t_start, t_avg, x, y, min_dt);
            draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x02), t_avg, t_end, x, y, min_dt);
        }
    } else {
        // If neither point is in-bounds, subdivide and iterate
        double t_avg = (t_end + t_start)/2;
        eval_func(t_avg, &x, &y, expr, stackpos);
        if (hypot((x - xi)*xscale, (y - yi)*yscale) < PLOT_MAX_ARC_LENGTH) return;
        //printf("Neither point in bounds on interval %f, %f -> (%f, %f), (%f, %f), %d\n", t_start, t_end, x, y, xi, yi, IN_BOUNDS(x, y));
        if (IN_BOUNDS(x, y)) {
            // [t_start, t_avg]
            draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x01) | 0x02, t_start, t_avg, x, y, min_dt);
            // [t_avg, t_end]
            draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x02) | 0x01, t_avg, t_end, x, y, min_dt);
        } else {
            draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x01), t_start, t_avg, x, y, min_dt);
            draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x02), t_avg, t_end, x, y, min_dt);
        }
    }
}

void draw_function_constant_ds(expression *expr, file_data *fd, uint8_t *color, cairo_t *cr) {
    SET_COLOR(cr, color);
    double x, y, xp, yp, dt, t, t_end;
    double *stackpos = fd->stack + fd->n_stack;
    t = 0;
    t_end = 1;
    double min_dt = PLOT_MIN_DT;
    uint32_t type = eval_func(0, &xp, &yp, expr, stackpos);
    if ((type & TYPE_MASK) != TYPE_POINT) {
        t = window_x0;
        t_end = window_x1;
        min_dt = PLOT_MAX_ARC_LENGTH/xscale;
    }
    nfev = 0;
    niev = 0;
    uint8_t flags;
    uint32_t n = (type>>8)/GET_STEP(type);
    //printf("Number of functions to be plotted: %d\n", n);
    expr->cache_size = 0;
    for (uint32_t i=0; i < n; i++) {
        flags = 0;
        eval_index = i;
        eval_func(t, &xp, &yp, expr, stackpos);
        if (IN_BOUNDS(xp, yp)) flags |= 0x01;
        eval_func(t_end, &x, &y, expr, stackpos);
        if (IN_BOUNDS(x, y)) flags |= 0x02;
        //printf("Initial values: (%f, %f) and (%f, %f), flags %02x\n", x, y, xp, yp, flags);
        draw_function_constant_ds_rec(expr, fd, cr, flags, t, t_end, xp, yp, min_dt);
    }
}

void draw_implicit_rec(expression *expr, file_data *fd, cairo_t *cr, uint8_t *color, double *area, int divisions) {
    /*fd->variable_list[0].pointer = area;
    fd->variable_list[0].type = 1<<8;
    fd->variable_list[1].pointer = area+2;
    fd->variable_list[1].type = 1<<8;*/
    double *stackpos = fd->stack + fd->n_stack;
    double *lstackpos = fd->lstack;
    eval_inter_2d(area, area+2, expr, stackpos, lstackpos);
    //expr->func->inter(expr->func, stackpos, lstackpos);
    //niev++;
    function *func = expr->func;
#ifdef PLOT_USE_INEQUALITY
    uint8_t ineq = 0;
    uint8_t singular = 0;
    if (func->oper == func_compare_sub_single) ineq = 1;
    if (func->oper == func_compare_sub) ineq = 1;
#endif
    if (lstackpos[0] > 0) {
        // Everything in the interval is greater than zero
#ifdef PLOT_USE_INEQUALITY
        if (ineq) {
            SET_COLOR_OPACITY(cr, color, 0.5);
            cairo_rectangle(cr, SCALE_XK(area[0]), SCALE_YK(area[3]), xscale*(area[1] - area[0]), yscale*(area[3] - area[2]));
            cairo_fill(cr);
            SET_COLOR_OPACITY(cr, color, 1);
        }
#endif
        return;
    }
    if (stackpos[0] < 0) {
        // No contours can exist in the given area
        //printf("    No contours found in ([%f, %f], [%f, %f])\n", area[0], area[1], area[2], area[3]);
        //cairo_rectangle(cr, SCALE_XK(area[0]), SCALE_YK(area[3]), xscale*(area[1] - area[0]), yscale*(area[3] - area[2]));
        //cairo_stroke(cr);
        return;
    }
    if (divisions >= PLOT_IMPLICIT_MAXDEPTH) {
        // Maximum depth reached
        //printf("    maximum depth reached on ([%f, %f], [%f, %f]) --> [%f, %f]\n", area[0], area[1], area[2], area[3], lstackpos[0], stackpos[0]);
        //cairo_rectangle(cr, SCALE_XK(area[0]), SCALE_YK(area[3]), xscale*(area[1] - area[0]), yscale*(area[3] - area[2]));
        //cairo_stroke(cr);
        function *arg1 = expr->func->first_arg;
        function *arg2 = arg1->next_arg;
        double x0 = area[0], x1 = area[1], y0 = area[2], y1 = area[3];
        double e00, e10, e11, e01;
        eval_func_2d(&x0, &y0, expr, stackpos); e00 = stackpos[0];
        eval_func_2d(&x1, &y0, expr, stackpos); e10 = stackpos[0];
        eval_func_2d(&x1, &y1, expr, stackpos); e11 = stackpos[0];
        eval_func_2d(&x0, &y1, expr, stackpos); e01 = stackpos[0];

        if ((e00 > 0) && (e10 > 0) && (e11 > 0) && (e01 > 0)) {
            // Interval calculation was wrong
#ifdef PLOT_USE_INEQUALITY
            if (ineq) {
                SET_COLOR_OPACITY(cr, color, 0.5);
                cairo_rectangle(cr, SCALE_XK(area[0]), SCALE_YK(area[3]), xscale*(area[1] - area[0]), yscale*(area[3] - area[2]));
                cairo_fill(cr);
                SET_COLOR_OPACITY(cr, color, 1);
            }
#endif
            return;
        }
        if (isinf(stackpos[0]) || isinf(lstackpos[0])) {
            // Singularity in the interval
            //printf("Singularity detected in ([%f, %f], [%f, %f])\n", area[0], area[1], area[2], area[3]);
            e00 = 1/e00;
            e10 = 1/e10;
            e11 = 1/e11;
            e01 = 1/e01;
            singular = 1;
        }
        double stroke_opacity = 1-singular;
        SET_COLOR_OPACITY(cr, color, stroke_opacity);

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
        if (((singular) || (edges == 0)) && (divisions < PLOT_IMPLICIT_HARDMAX)) {
            //printf("Unresolved detail in ([%f, %f], [%f, %f])\n", area[0], area[1], area[2], area[3]);
            // Subdivide
            double x0 = area[0], x1 = area[1], y0 = area[2], y1 = area[3];
            double xm = (x0 + x1)/2, ym = (y0 + y1)/2;
            //printf("subdividing ([%f, %f], [%f, %f]), %f, %f\n", x0, x1, y0, y1, xm, ym);
            area[1] = xm; area[3] = ym;
            draw_implicit_rec(expr, fd, cr, color, area, divisions+1);
            area[0] = xm; area[1] = x1; area[2] = y0; area[3] = ym;
            draw_implicit_rec(expr, fd, cr, color, area, divisions+1);
            area[0] = x0; area[1] = xm; area[2] = ym; area[3] = y1;
            draw_implicit_rec(expr, fd, cr, color, area, divisions+1);
            area[0] = xm; area[1] = x1; area[2] = ym; area[3] = y1;
            draw_implicit_rec(expr, fd, cr, color, area, divisions+1);
            return;
        }
        // If two edges are selected, then draw the line between those two edges. If three 
        // edges are selected, then one of the corners must be zero so the contour must go
        // through one of the corners. If four edges are selected, there are two contour
        // lines, each going through opposite sides of the area.
        double sx0 = SCALE_XK(x0), sx1 = SCALE_XK(x1), sy0 = SCALE_YK(y0), sy1 = SCALE_YK(y1);
        double snpos = SCALE_XK(npos), sepos = SCALE_YK(epos), sspos = SCALE_XK(spos), swpos = SCALE_YK(wpos);
        // If the contout passes through two corners, only draw one line
        if ((edges == 0xf) && (((e00 == 0) && (e11 == 0)) || ((e10 == 0) && (e01 == 0)))) edges = 0x5;
        switch (edges) {
            case 0xf:
                // If all four edges are selected, we know that all of the corners have
                // different values from their neighbors. No more than two corners can
                // be zero and they must be on opposite corners. If we sample two adjacent
                // corners, at least one will not be zero
#ifdef PLOT_USE_INEQUALITY
                if (ineq) {
                    SET_COLOR_OPACITY(cr, color, 0.5);
                    if ((e00 > 0) || (e10 < 0)) {
                        cairo_move_to(cr, sx0, sy0);
                        cairo_line_to(cr, sspos, sy0);
                        cairo_line_to(cr, snpos, sy1);
                        cairo_line_to(cr, sx1, sy1);
                        cairo_line_to(cr, sx1, sepos);
                        cairo_line_to(cr, sx0, swpos);
                    } else {
                        cairo_move_to(cr, sx1, sy0);
                        cairo_line_to(cr, sspos, sy0);
                        cairo_line_to(cr, snpos, sy1);
                        cairo_line_to(cr, sx0, sy1);
                        cairo_line_to(cr, sx0, sepos);
                        cairo_line_to(cr, sx1, swpos);
                    }
                    cairo_close_path(cr);
                    cairo_fill(cr);
                }
                SET_COLOR_OPACITY(cr, color, stroke_opacity);
#endif
                cairo_move_to(cr, sspos, sy0);
                cairo_line_to(cr, snpos, sy1);
                cairo_stroke(cr);
                cairo_move_to(cr, sx0, swpos);
                cairo_line_to(cr, sx1, sepos);
                cairo_stroke(cr);
                break;
            case 0x7:
            case 0xd:
            case 0x5:
                cairo_move_to(cr, sspos, sy0);
                cairo_line_to(cr, snpos, sy1);
#ifdef PLOT_USE_INEQUALITY
                if (ineq) {
                    cairo_stroke_preserve(cr);
                    SET_COLOR_OPACITY(cr, color, 0.5);
                    if ((e00 > 0) || (e10 < 0)) {
                        cairo_line_to(cr, sx0, sy1);
                        cairo_line_to(cr, sx0, sy0);
                    } else {
                        cairo_line_to(cr, sx1, sy1);
                        cairo_line_to(cr, sx1, sy0);
                    }
                    cairo_close_path(cr);
                    cairo_fill(cr);
                    SET_COLOR_OPACITY(cr, color, stroke_opacity);
                } else
#endif
                cairo_stroke(cr);
                break;
            case 0xb:
            case 0xe:
            case 0xa:
                cairo_move_to(cr, sx0, swpos);
                cairo_line_to(cr, sx1, sepos);
#ifdef PLOT_USE_INEQUALITY
                if (ineq) {
                    cairo_stroke_preserve(cr);
                    SET_COLOR_OPACITY(cr, color, 0.5);
                    if ((e00 > 0) || (e01 < 0)) {
                        cairo_line_to(cr, sx1, sy0);
                        cairo_line_to(cr, sx0, sy0);
                    } else {
                        cairo_line_to(cr, sx1, sy1);
                        cairo_line_to(cr, sx0, sy1);
                    }
                    cairo_close_path(cr);
                    cairo_fill(cr);
                    SET_COLOR_OPACITY(cr, color, stroke_opacity);
                } else
#endif
                cairo_stroke(cr);
                break;
            case 0x3:
                cairo_move_to(cr, sx1, sepos);
                cairo_line_to(cr, snpos, sy1);
#ifdef PLOT_USE_INEQUALITY
                if (ineq) {
                    cairo_stroke_preserve(cr);
                    SET_COLOR_OPACITY(cr, color, 0.5);
                    if ((e11 > 0) || (e01 < 0)) {
                        cairo_line_to(cr, sx1, sy1);
                    } else {
                        cairo_line_to(cr, sx0, sy1);
                        cairo_line_to(cr, sx0, sy0);
                        cairo_line_to(cr, sx1, sy0);
                    }
                    cairo_close_path(cr);
                    cairo_fill(cr);
                    SET_COLOR_OPACITY(cr, color, stroke_opacity);
                } else
#endif
                cairo_stroke(cr);
                break;
            case 0x6:
                cairo_move_to(cr, sx1, sepos);
                cairo_line_to(cr, sspos, sy0);
#ifdef PLOT_USE_INEQUALITY
                if (ineq) {
                    cairo_stroke_preserve(cr);
                    SET_COLOR_OPACITY(cr, color, 0.5);
                    if ((e10 > 0) || (e00 < 0)) {
                        cairo_line_to(cr, sx1, sy0);
                    } else {
                        cairo_line_to(cr, sx0, sy0);
                        cairo_line_to(cr, sx0, sy1);
                        cairo_line_to(cr, sx1, sy1);
                    }
                    cairo_close_path(cr);
                    cairo_fill(cr);
                    SET_COLOR_OPACITY(cr, color, stroke_opacity);
                } else
#endif
                cairo_stroke(cr);
                break;
            case 0xc:
                cairo_move_to(cr, sx0, swpos);
                cairo_line_to(cr, sspos, sy0);
#ifdef PLOT_USE_INEQUALITY
                if (ineq) {
                    cairo_stroke_preserve(cr);
                    SET_COLOR_OPACITY(cr, color, 0.5);
                    if ((e10 < 0) || (e00 > 0)) {
                        cairo_line_to(cr, sx0, sy0);
                    } else {
                        cairo_line_to(cr, sx1, sy0);
                        cairo_line_to(cr, sx1, sy1);
                        cairo_line_to(cr, sx0, sy1);
                    }
                    cairo_close_path(cr);
                    cairo_fill(cr);
                    SET_COLOR_OPACITY(cr, color, stroke_opacity);
                } else
#endif
                cairo_stroke(cr);
                break;
            case 0x9:
                cairo_move_to(cr, sx0, swpos);
                cairo_line_to(cr, snpos, sy1);
#ifdef PLOT_USE_INEQUALITY
                if (ineq) {
                    cairo_stroke_preserve(cr);
                    SET_COLOR_OPACITY(cr, color, 0.5);
                    if ((e01 > 0) || (e00 < 0)) {
                        cairo_line_to(cr, sx0, sy1);
                    } else {
                        cairo_line_to(cr, sx1, sy1);
                        cairo_line_to(cr, sx1, sy0);
                        cairo_line_to(cr, sx0, sy0);
                    }
                    cairo_close_path(cr);
                    cairo_fill(cr);
                    SET_COLOR_OPACITY(cr, color, stroke_opacity);
                } else
#endif
                cairo_stroke(cr);
                break;
            default:
                break;
        }
    } else {
        // Subdivide
        double x0 = area[0], x1 = area[1], y0 = area[2], y1 = area[3];
        double xm = (x0 + x1)/2, ym = (y0 + y1)/2;
        //printf("subdividing ([%f, %f], [%f, %f]), %f, %f\n", x0, x1, y0, y1, xm, ym);
        area[1] = xm; area[3] = ym;
        draw_implicit_rec(expr, fd, cr, color, area, divisions+1);
        area[0] = xm; area[1] = x1; area[2] = y0; area[3] = ym;
        draw_implicit_rec(expr, fd, cr, color, area, divisions+1);
        area[0] = x0; area[1] = xm; area[2] = ym; area[3] = y1;
        draw_implicit_rec(expr, fd, cr, color, area, divisions+1);
        area[0] = xm; area[1] = x1; area[2] = ym; area[3] = y1;
        draw_implicit_rec(expr, fd, cr, color, area, divisions+1);
    }
}

void draw_implicit(expression *expr, file_data *fd, uint8_t *color, cairo_t *cr) {
    SET_COLOR(cr, color);
    double temp[4] = {window_x0, window_x1, window_y0, window_y1};
    uint32_t (*old_oper)(void*, double*) = expr->func->oper;
    if (expr->func->oper == func_equals) expr->func->oper = func_sub;
    else if (expr->func->oper == func_compare_single) expr->func->oper = func_compare_sub_single;
    else if (expr->func->oper == func_compare) expr->func->oper = func_compare_sub;
    draw_implicit_rec(expr, fd, cr, color, temp, 0);
    expr->func->oper = old_oper;
}

void reset_overlay(file_data *fd) {
    cairo_t *cr = cairo_create(fd->overlay);
    cairo_operator_t oper = cairo_get_operator(cr);
    cairo_set_operator(cr, CAIRO_OPERATOR_CLEAR);
    cairo_set_source_rgba(cr, 0, 0, 0, 0);
    cairo_paint(cr);
    cairo_set_operator(cr, oper);
    fd->use_overlay = 0;
    fd->click_expr = NULL;
    cairo_destroy(cr);
}

uint8_t find_nearest_point(file_data *fd, expression *expr, double x, double y, double *mindist, int *index, double *xf, double *yf, uint8_t force) {
    double *stack = fd->stack, *lstack = fd->lstack;
    uint32_t n_stack = fd->n_stack;
    double dist;
    *mindist = INFINITY;
    if (expr->flags & EXPRESSION_PLOTTABLE) {
        if ((expr->flags & EXPRESSION_FIXED) && ((expr->value_type & TYPE_MASK) == TYPE_POINT)) {
            int len = (expr->value_type) >> 8;
            double *ptr = expr->value;
            double pt_x, pt_y;
            for (int p=0; p < len; p+=2) {
                pt_x = SCALE_XK(ptr[p]);
                pt_y = SCALE_YK(ptr[p+1]);
                dist = hypot(pt_x - x, pt_y - y);
                if ((p == 0) || (dist < *mindist)) {
                    *mindist = dist;
                    *index = p/2;
                    *xf = ptr[p];
                    *yf = ptr[p+1];
                }
            }
        } else if ((expr->func->oper == func_equals) || (expr->func->oper == func_compare) || (expr->func->oper == func_compare_single)) {
            double x0 = x/xscale + window_x0;
            double y0 = window_y1 - y/yscale;
            double xc = x0, yc = y0;
            double temp, f00, f10, f01;
            double inter[4] = {x0, 0, y0, 0};
            inter[1] = inter[0] + CLICK_MAX_RADIUS/xscale;
            inter[0] -= CLICK_MAX_RADIUS/xscale;
            inter[3] = inter[2] + CLICK_MAX_RADIUS/yscale;
            inter[2] -= CLICK_MAX_RADIUS/yscale;
            eval_inter_2d(inter, inter+2, expr, stack+n_stack, lstack);
            uint32_t (*old_oper)(void*, double*) = expr->func->oper;
            if (expr->func->oper == func_equals) expr->func->oper = func_sub;
            else if (expr->func->oper == func_compare_single) expr->func->oper = func_compare_sub_single;
            else if (expr->func->oper == func_compare) expr->func->oper = func_compare_sub;
            if (((stack[n_stack] >= 0) && (lstack[0] <= 0)) || force) {
                for (int j=0; j < CLICK_MAX_ITER; j++) {
                    eval_func_2d(&x0, &y0, expr, stack+n_stack);
                    f00 = stack[n_stack];
                    temp = x0+PLOT_DT_DERIV;
                    eval_func_2d(&temp, &y0, expr, stack+n_stack);
                    f10 = stack[n_stack];
                    temp = y0+PLOT_DT_DERIV;
                    eval_func_2d(&x0, &temp, expr, stack+n_stack);
                    f01 = stack[n_stack];
                    temp = f00 / (pow((f10 - f00)/PLOT_DT_DERIV, 2) + pow((f01 - f00)/PLOT_DT_DERIV, 2));
                    x0 -= temp*(f10 - f00)/PLOT_DT_DERIV;
                    y0 -= temp*(f01 - f00)/PLOT_DT_DERIV;
                }
                *mindist = hypot((x0-xc)*xscale, (y0-yc)*yscale);
                *index = 0;
                *xf = x0;
                *yf = y0;
                eval_func_2d(&x0, &y0, expr, stack+n_stack);
                if (fabs(stack[n_stack]) >= CLICK_IMPLICIT_MAX_ERROR) {
                    expr->func->oper = old_oper;
                    return FIND_NEAREST_NOCONV;
                }
            }
            expr->func->oper = old_oper;
        } else if (!(expr->flags & EXPRESSION_FIXED) && (expr->cache_size != -1)) {
            for (uint32_t i=0; i < expr->cache_size; i+=2) {
                dist = hypot(SCALE_XK(expr->value[i]) - x, SCALE_YK(expr->value[i+1]) - y);
                if (dist < *mindist) {
                    *mindist = dist;
                    *xf = expr->value[i];
                    *yf = expr->value[i+1];
                }
            }
        } else {
            return FIND_NEAREST_UNSUPPORTED;
        }
    }
    return 0;
}

void draw_labeled_point(cairo_t *cr, double xf, double yf) {
    cairo_set_source_rgb(cr, 0, 0, 1);
    cairo_set_line_width(cr, 1);
    double xk = SCALE_XK(xf), yk = SCALE_YK(yf);
    cairo_rectangle(cr, xk - CLICK_MAX_RADIUS, yk - CLICK_MAX_RADIUS, 2*CLICK_MAX_RADIUS, 2*CLICK_MAX_RADIUS);
    cairo_stroke(cr);
    
    cairo_select_font_face(cr, "Verdana", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    cairo_move_to(cr, xk, yk);
    char temp[50];
    sprintf(temp, "(%f, %f)", xf, yf);
    cairo_show_text(cr, temp);
}

static gboolean button_press_callback (GtkWidget *event_box, GdkEventButton *event, gpointer data_pointer) {
    file_data *fd = (file_data*)(data_pointer);
    g_print ("Clicked at %f, %f\n", event->x, event->y);
    expression *expression_list = fd->expression_list, *expr;
    double dist, mindist=INFINITY;
    int expr_clicked=0, index=0;
    int temp_index;
    double *stack = fd->stack, *lstack = fd->lstack;
    uint32_t n_stack = fd->n_stack;
    double xf, yf, temp_xf, temp_yf;
    for (int i=0; i < fd->n_expr; i++) {
        expr = expression_list+i;
        find_nearest_point(fd, expr, event->x, event->y, &dist, &temp_index, &temp_xf, &temp_yf, 0);
        if ((i == 0) || (dist < mindist)) {
            xf = temp_xf;
            yf = temp_yf;
            index = temp_index;
            mindist = dist;
            expr_clicked = i;
        }
    }
    if ((mindist != INFINITY) && (mindist < CLICK_MAX_RADIUS)) {
        printf("User clicked on expression %d, index %d\n", expr_clicked, index);
        cairo_t *cr = cairo_create(fd->overlay);
        draw_labeled_point(cr, xf, yf);
        cairo_destroy(cr);
        fd->use_overlay = 1;
        fd->click_expr = expression_list + expr_clicked;
        gtk_widget_queue_draw(fd->drawing_area);
    }
    click_x = event->x;
    click_y = event->y;
    click_state = 1;
    return TRUE;
}

static gboolean button_release_callback (GtkWidget *event_box, GdkEventButton *event, gpointer data_pointer) {
    file_data *fd = (file_data*)(data_pointer);
    printf("Clearing overlay\n");
    reset_overlay(fd);
    click_x = event->x;
    click_y = event->y;
    click_state = 0;
    gtk_widget_queue_draw(fd->drawing_area);
    return TRUE;
}

void redraw_all(file_data *fd) {
    reset_overlay(fd);

    cairo_t *cr = cairo_create(fd->surface);
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_paint(cr);
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
    double x0_scaled = ticksize*floor(window_x0/ticksize);
    for (double tick=x0_scaled; tick <= window_x1; tick+=ticksize) {
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
    double y0_scaled = ticksize*floor(window_y0/ticksize);
    for (double tick=y0_scaled; tick <= window_y1; tick+=ticksize) {
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
    uint32_t color_type, color_len;
    for (int i=0; i < fd->n_expr; i++) {
        expr = expression_list+i;
        if (expr->color_pointer) {
            printf("Expression %p (%d) has color pointer\n", expr, i+1);
            color_type = expr->color_pointer->oper(expr->color_pointer, stack+n_stack);
            if (!((color_type & TYPE_MASK) == TYPE_COLOR)) {
                printf("ERROR: Color expression must return color\n");
                exit(EXIT_FAILURE);
            }
            print_object(color_type, stack+n_stack); printf("\n");
            color_len = color_type>>8;
            fd->n_stack += color_len;
            n_stack += color_len;
        } else color_len = 0;
        if (expr->flags & EXPRESSION_PLOTTABLE) {
            if ((expr->flags & EXPRESSION_FIXED) && ((expr->value_type & TYPE_MASK) == TYPE_POINT)) {
#ifdef DEBUG_PLOT
                t3 = clock();
#endif
                int len = (expr->value_type) >> 8;
                double *ptr = expr->value;
                SET_COLOR(cr, expr->color);
                for (int p=0; p < len; p+=2) {
                    if (expr->color_pointer) SET_COLOR(cr, (stack+n_stack-color_len+(p/2*3)%color_len));
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
                    //if (expr->color_pointer) SET_COLOR(cr, (expression_list[20].value+color_pos));
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
            } else if ((expr->func->oper == func_equals) || (expr->func->oper == func_compare) || (expr->func->oper == func_compare_single)) {
                if (!(expr->func->inter)) {
                    printf("ERROR: interval function needed for implicit plotting\n");
                    continue;
                }
#ifdef DEBUG_PLOT
                t3 = clock();
                nfev = 0;
                niev = 0;
#endif
                draw_implicit(expression_list+i, fd, expr->color, cr);
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
                draw_function_constant_ds(expression_list+i, fd, expr->color, cr);
#ifdef DEBUG_PLOT
                t4 = clock();
                printf("Plotted expression %p (%d) in %luus, nfev: %d, niev: %d, interval function %p\n", expression_list+i, i+1, t4-t3, nfev, niev, expression_list[i].func->inter);
#endif
            }
        }
        n_stack -= color_len;
        fd->n_stack -= color_len;
    }
    cairo_destroy(cr);
    t2 = clock();
    g_print("Redraw took %luus, %d expressions, bounds %f %f %f %f\n", t2-t1, n_expr, window_x0, window_y0, window_x1, window_y1);
}

gboolean redraw_callback(GtkWidget *widget, cairo_t *cr, gpointer data_pointer) {
    file_data *fd = (file_data*)(data_pointer);
    cairo_set_source_surface(cr, fd->surface, 0, 0);
    cairo_paint(cr);
    if (fd->use_overlay) {
        cairo_set_source_surface(cr, fd->overlay, 0, 0);
        cairo_paint(cr);
    }
    return FALSE;
}

void run_action(file_data *fd, function *action) {
    clock_t t3 = clock();
    action->oper(action, (fd->stack)+(fd->n_stack));
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
        redraw_all(fd);
        gtk_widget_queue_draw(fd->drawing_area);
    }
}

static gboolean configure_callback(GtkWidget *widget, GdkEventConfigure *event, gpointer data_pointer) {
    file_data *fd = (file_data*)(data_pointer);
    if (fd->surface) cairo_surface_destroy(fd->surface);
    fd->surface = gdk_window_create_similar_surface(gtk_widget_get_window(widget), CAIRO_CONTENT_COLOR,
                                                    gtk_widget_get_allocated_width(widget), gtk_widget_get_allocated_height(widget));
    if (fd->overlay) cairo_surface_destroy(fd->overlay);
    fd->overlay = gdk_window_create_similar_surface(gtk_widget_get_window(widget), CAIRO_CONTENT_COLOR_ALPHA,
                                                    gtk_widget_get_allocated_width(widget), gtk_widget_get_allocated_height(widget));
    redraw_all(fd);
    return TRUE;
}

static gboolean scroll_callback (GtkWidget *event_box, GdkEventScroll *event, gpointer data_pointer) {
    file_data *fd = (file_data*)(data_pointer);
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
        redraw_all(fd);
        gtk_widget_queue_draw(fd->drawing_area);
    }
    return TRUE;
}

static gboolean motion_callback(GtkWidget *event_box, GdkEventMotion *event, gpointer data_pointer) {
    file_data *fd = (file_data*)(data_pointer);
    if ((click_state) && (fd->use_overlay) && (fd->click_expr)) {
        double xf, yf, dist;
        int index;
        uint8_t status = find_nearest_point(fd, fd->click_expr, event->x, event->y, &dist, &index, &xf, &yf, 1);
        cairo_t *cr = cairo_create(fd->overlay);
        // Clear background
        cairo_operator_t oper = cairo_get_operator(cr);
        cairo_set_operator(cr, CAIRO_OPERATOR_CLEAR);
        cairo_set_source_rgba(cr, 0, 0, 0, 0);
        cairo_paint(cr);
        cairo_set_operator(cr, oper);
        // Draw rectangle
        if (!status) draw_labeled_point(cr, xf, yf);
        cairo_destroy(cr);
        gtk_widget_queue_draw(fd->drawing_area);
    } else if (click_state) {
        double dx = (event->x - click_x)/xscale;
        double dy = (click_y - event->y)/yscale;
        window_x0 -= dx; window_x1 -= dx; window_y0 -= dy; window_y1 -= dy;
        click_x = event->x;
        click_y = event->y;
        uint8_t color[4] = {255, 255, 255, 0};
        redraw_all(fd);
        gtk_widget_queue_draw(fd->drawing_area);
    }
    return TRUE;
}

static gboolean timeout_callback(gpointer data_pointer) {
    printf("timeout_callback\n");
    clock_t t3 = clock();
    file_data *fd = (file_data*)(data_pointer);
    run_action(fd, fd->expression_list[ticker_target].func);
    if (run_ticker) return TRUE;
    return FALSE;
}

static gboolean keypress_callback(GtkWidget *widget, GdkEventKey *event, gpointer data_pointer) {
    file_data *fd = (file_data*)data_pointer;
    printf("Event is %d\n", event->keyval);
    if (event->keyval == 's') {
        run_ticker = 1;
        if (ticker_target >= 0) g_timeout_add(ticker_step, timeout_callback, data_pointer);
    } else if (event->keyval == 'e') {
        run_ticker = 0;
    } else if (event->keyval == 'i') {
        treeview_activate(fd);
        printf("inspecting\n");
    }
    expression *expr;
    function *action;
    for (int i=0; i < fd->n_expr; i++) {
        expr = fd->expression_list+i;
        if ((expr->func) && (expr->func->oper == func_onkeypress) && (((double*)(expr->value))[0] == event->keyval)) {
            action = *((function**)(((double*)(expr->value))+1));
            printf("Calling action %p\n", action);
            run_action(fd, action);
        }
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
    g_signal_connect(G_OBJECT(drawing_area), "draw", G_CALLBACK(redraw_callback), fd);
    g_signal_connect(G_OBJECT(drawing_area), "configure-event", G_CALLBACK(configure_callback), fd);
    g_signal_connect(window, "key-press-event", G_CALLBACK(keypress_callback), fd);
    gtk_container_add(GTK_CONTAINER(event_box), drawing_area);

    gtk_widget_show_all(window);
}

int main (int argc, char **argv) {
    memset(variable_list, 0, 10*sizeof(variable));
    memset(expression_list, 0, 10*sizeof(expression));
    printf("First stack positions %p, %p\n", stack, stack+1);
    variable_list[0] = new_variable("x", 0, VARIABLE_IN_SCOPE | VARIABLE_INTERVAL, NULL);
    variable_list[1] = new_variable("y", 0, VARIABLE_IN_SCOPE | VARIABLE_INTERVAL, NULL);
    double pi = M_PI;
    variable_list[2] = new_variable("\\pi", 1<<8, VARIABLE_IN_SCOPE, &pi);
    variable_list[3] = new_variable("\\theta", 0, VARIABLE_IN_SCOPE, NULL);

    fd.expression_list = expression_list;
    fd.variable_list = variable_list;
    fd.function_list = function_list;
    fd.stack = stack;
    fd.lstack = lstack;
    fd.n_expr = 0;
    fd.n_var = 4;
    fd.n_func = 0;
    fd.n_stack = 0;
    fd.deptable = deptable;
    fd.surface = NULL;
    fd.overlay = NULL;
    fd.use_overlay = 0;
    fd.click_expr = NULL;
    
    uint32_t n_func = 0;
    uint32_t n_var = 0;
    load_file(argv[1], &fd);
    top_expr = parse_file(&fd, stringbuf);
    printf("Parsing completed. %d function blocks, %d variables, %d expressions\n", fd.n_func, fd.n_var, fd.n_expr);
    ticker_target = (int)parse_double(argv[2]);
    ticker_step = (int)parse_double(argv[3]);
    printf("Expression %d will be evaluated every %d milliseconds\n", ticker_target, ticker_step);
    //double area[4] = {0.9, 1, 1.2, 1.25};
    //eval_inter_2d(area, area+2, expression_list+5, stack+n_stack, lstack);
    //printf("interval is [%f, %f]\n", lstack[0], stack[n_stack]);
    //exit(EXIT_FAILURE);

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
