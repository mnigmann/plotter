#include <gtk/gtk.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#include "parse.h"
#include "functions.h"
#include "linalg_functions.h"
#include "config.h"
#include "treeview.h"

#define SCALE_XK(v) ((v - window_x0)*xscale)
#define SCALE_YK(v) ((window_y1 - v)*yscale)

#define CLIP_WIDTH(v) (v<0 ? 0 : (v>WIDTH ? WIDTH : v))
#define CLIP_HEIGHT(v) (v<0 ? 0 : (v>HEIGHT ? HEIGHT : v))
#define IN_BOUNDS(x, y) ((window_x0 <= x) && (x <= window_x1) && (window_y0 <= y) && (y <= window_y1) && (x == x) && (y == y))
#define SET_COLOR(cr, color) cairo_set_source_rgb(cr, color[0]/256.0, color[1]/256.0, color[2]/256.0)
#define SET_COLOR_HEX(cr, color) cairo_set_source_rgb(cr, ((color>>16)&255)/256.0, ((color>>8)&255)/256.0, (color&255)/256.0)
#define SET_COLOR_OPACITY(cr, color, opacity) cairo_set_source_rgba(cr, color[0]/256.0, color[1]/256.0, color[2]/256.0, opacity)

#define GET_STEP(type) (step_table[(type) & TYPE_MASK])
const static uint32_t step_table[8] = {1, 2, 3, 1, 1, MAX_POLYGON_SIZE*2, 0, 0};

double window_x0 = -10;
double window_y0 = -10;
double window_x1 = 10;
double window_y1 = 10;
double xscale;
double yscale;
uint16_t nfev = 0;
uint16_t niev = 0;
double enclosed_area = 0;

function function_list[2048];
variable variable_list[512];
expression expression_list[100];
expression *deptable[300];
double stack[65536];
double lstack[65536];
char stringbuf[1024];
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

int da_width = INIT_WIDTH, da_height = INIT_HEIGHT;
#define WIDTH da_width
#define HEIGHT da_height

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

void find_extremum(expression *expr, double *stackpos, double x, double y, double xp, double yp, double xpp, double ypp, uint8_t mode) {
    double u1, v1, u2, v2, ex, temp, val, m1, m2;
    double y_ex, ydiff;
    for (int i=0; i < 10; i++) {
        u1 = x - xp; v1 = y - yp;
        u2 = xpp - xp; v2 = ypp - yp;
        temp = fabs(yp) * PLOT_EXTREMUM_EPS;
        if ((fabs(v1) < temp) && (fabs(v2) < temp)) break;
        ex = (u2*u2*v1 - u1*u1*v2)/(2*(u2*v1 - u1*v2));
        eval_func(xp+ex, &temp, &y_ex, expr, stackpos);
        if (mode) ydiff = y_ex - yp;
        else ydiff = yp - y_ex;
        // We know without loss of generality that y > yp < ypp
        // If ex > 0 then x_ex is between xp and x
        //     If y_ex is between yp and y, then we use it as an endpoint
        //     If y_ex is less than than yp, then xp = x_ex and xpp = x_p
        // If ex < 0 then x_ex is between xp and xpp
        if (ex == 0) {
            x = (x + xp)/2;
            eval_func(x, &temp, &y, expr, stackpos);
            xpp = (xpp + xp)/2;
            eval_func(xpp, &temp, &ypp, expr, stackpos);
        } else if ((ex > 0) == (xpp > xp)) {
            // towards xpp
            if (ydiff > 0) {
                // New point is more extreme than the current extremum
                x = xp;
                y = yp;
                xp += ex;
                yp = y_ex;
            } else {
                // New point is less extreme than the current extremum
                xpp = xp + ex;
                ypp = y_ex;
            }
        } else {
            // towards x
            if (ydiff > 0) {
                xpp = xp;
                ypp = yp;
                xp += ex;
                yp = y_ex;
            } else {
                x = xp + ex;
                y = y_ex;
            }
        }
        printf("    iterated extremum at %f, value %f, ex %f; (%f, %f) (%f, %f) (%f, %f)\n", xp, yp, ex, xpp, ypp, xp, yp, x, y);
    }
    uint32_t old_n = expr->n_special_points;
    expr->special_points = realloc(expr->special_points, 3*sizeof(double)*(old_n + 1));
    expr->special_points[3*old_n] = xp;
    expr->special_points[3*old_n+1] = xp;
    expr->special_points[3*old_n+2] = yp;
    expr->n_special_points++;
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
        //printf("        iterated discontinuity %f -> (%f, %f) and %f -> (%f, %f)\n", tp, xp, yp, t, x, y);
    }
    //printf("    function is (%f, %f)\n", xp, yp);
    //printf("    broken to (%f, %f)\n", x, y);
    cairo_line_to(cr, SCALE_XK(xp), SCALE_YK(yp));
    cairo_move_to(cr, SCALE_XK(x), SCALE_YK(y));
}

void draw_function_constant_ds_rec(expression *expr, file_data *fd, cairo_t *cr, uint8_t flags, double t_start, double t_end, double xi, double yi, double min_dt) {
    double x, y, xp, yp, xpp, ypp, tp, tpp, dt, ds, dsp;
    //printf("t_end %f, t_start %f, min_dt %f, flags %02x, interval %p\n", t_end, t_start, min_dt, flags, expr->func->inter);
    double *stackpos = fd->stack + fd->n_stack;
    double *lstackpos = fd->lstack;
    if (t_end - t_start < min_dt) {
        // Final check if there could be another line to draw.
        // The two endpoints must both be off-screen, but not on the 
        // same side.
        if ((flags & 0x0f) && (flags & 0xf0) && !((flags>>4) & flags)) {
            //printf("Extra line needed between %f and %f, flags %02x\n", t_start, t_end, flags);
            eval_func(t_start, &xp, &yp, expr, stackpos);
            eval_func(t_end, &x, &y, expr, stackpos);
            // Truncate the line to fit in the screen. A line must cross
            // at least two distinct edges
            // Check t_start (xp, yp)
            double a;
            if (flags & 0x0c) {
                a = (flags & 0x08) ? window_y1 : window_y0;
                xp += (a - yp) * (x - xp) / (y - yp);
                yp = a;
                flags &= 0xfc;
                if (xp < window_x0) flags |= 0x01;
                if (xp > window_x1) flags |= 0x02;
            }
            // Check t_end (x, y)
            if (flags & 0xc0) {
                a = (flags & 0x80) ? window_y1 : window_y0;
                x += (a - y) * (x - xp) / (y - yp);
                y = a;
                flags &= 0xcf;
                if (xp < window_x0) flags |= 0x10;
                if (xp > window_x1) flags |= 0x20;
            }
            if (!((flags >> 4) & flags) && flags) {
                if (flags & 0x03) {
                    a = (flags & 0x01) ? window_x0 : window_x1;
                    yp += (a - xp) * (y - yp) / (x - xp);
                    xp = a;
                    flags &= 0xf0;
                }
                if (flags & 0x30) {
                    a = (flags & 0x10) ? window_x0 : window_x1;
                    y += (a - xp) * (y - yp) / (x - xp);
                    x = a;
                    flags &= 0x0f;
                }
            }
            cairo_move_to(cr, SCALE_XK(xp), SCALE_YK(yp));
            cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
            //printf("Line between (%f, %f) and (%f, %f)\n", xp, yp, x, y);
            cairo_stroke(cr);
        }
        //printf("Returning due to minimum dt\n");
        if ((flags & 0x0f) && (flags & 0xf0)) return;
    }
    //printf("recursive step from %f to %f, flags %02x\n", t_start, t_end, flags);
    if (!(flags & 0x0f)) {
        // If the lower point is in-bounds, start from there and go up
        eval_func(t_start, &xp, &yp, expr, stackpos);
        cairo_move_to(cr, SCALE_XK(xp), SCALE_YK(yp));
        eval_func(t_start+PLOT_DT_DERIV, &x, &y, expr, stackpos);
        //printf("intial point (%.12f, %.12f), (%.12f, %.12f), (%.12f, %.12f)\n", xpp, ypp, xp, yp, x, y);
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
            //printf("    function is (%f, %f) at %f\n", x, y, t_start);
            if ((x == x) && (y == y)) {
                ds = hypot((x-xp)*xscale, (y-yp)*yscale);
                dt = PLOT_MAX_ARC_LENGTH*dt/ds;
                if ((expr->flags & EXPRESSION_EXPLICIT) && (yp > ypp) && (yp > y) && (step >= 1)) {
                    printf("Possible maximum at %f -> (%f, %f)\n", tp, xp, yp);
                    find_extremum(expr, stackpos, t_start, y, tp, yp, tpp, ypp, 1);
                } else if ((expr->flags & EXPRESSION_EXPLICIT) && (yp < ypp) && (yp < y) && (step >= 1)) {
                    printf("Possible minimum at %f -> (%f, %f) > (%f, %f) < (%f, %f)\n", tp, x, y, xp, yp, xpp, ypp);
                    find_extremum(expr, stackpos, t_start, y, tp, yp, tpp, ypp, 0);
                }
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
                flags &= 0xf0;
                if (x < window_x0) flags |= 0x01;
                if (x > window_x1) flags |= 0x02;
                if (y < window_y0) flags |= 0x04;
                if (y > window_y1) flags |= 0x08;
                if (!((x == x) && (y == y))) flags |= 0x0f;
                cairo_stroke(cr);
                draw_function_constant_ds_rec(expr, fd, cr, flags, t_start, t_end, x, y, min_dt);
                break;
            }
        }
        /*eval_func(t_end, &x, &y, expr, stackpos);
        if ((x == x) && (y == y)) {
            cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
            step++;
        }*/
        cairo_stroke(cr);
    } else if (!(flags & 0xf0)) {
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
        uint32_t step = 0;
        while (t_end > t_start) {
            eval_func(t_end, &x, &y, expr, stackpos);
            //printf("    function is (%f, %f) at %f\n", x, y, t_start);
            if ((x == x) && (y == y)) {
                ds = hypot((x-xp)*xscale, (y-yp)*yscale);
                dt = PLOT_MAX_ARC_LENGTH*dt/ds;
                if ((expr->flags & EXPRESSION_EXPLICIT) && (yp > ypp) && (yp > y) && (step >= 1)) {
                    printf("Possible maximum at %f -> (%f, %f) < (%f, %f) > (%f, %f)\n", tp, x, y, xp, yp, xpp, ypp);
                    find_extremum(expr, stackpos, t_end, y, tp, yp, tpp, ypp, 1);
                } else if ((expr->flags & EXPRESSION_EXPLICIT) && (yp < ypp) && (yp < y) && (step >= 1)) {
                    printf("Possible minimum at %f -> (%f, %f) > (%f, %f) < (%f, %f)\n", tp, x, y, xp, yp, xpp, ypp);
                    find_extremum(expr, stackpos, t_end, y, tp, yp, tpp, ypp, 0);
                }
                if (ds > PLOT_MAX_ARC_LENGTH*PLOT_DISCONTINUITY_THRESHOLD) {
                    //printf("Possible discontinuity between %f and %f: (%f, %f) and (%f, %f)\n", t_end, tp, x, y, xp, yp);
                    find_discontinuity(expr, stackpos, t_end, x, y, tp, xp, yp, cr);
                }
                cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
                xpp = xp; ypp = yp; xp = x; yp = y; tpp = tp; tp = t_end; dsp = ds;
                t_end -= dt;
                step++;
            }
            if (!IN_BOUNDS(x, y)) {
                // Leaving the bounds, so iterate on [t, t_end]
                flags &= 0x0f;
                if (x < window_x0) flags |= 0x10;
                if (x > window_x1) flags |= 0x20;
                if (y < window_y0) flags |= 0x40;
                if (y > window_y1) flags |= 0x80;
                if (!((x == x) && (y == y))) flags |= 0xf0;
                cairo_stroke(cr);
                draw_function_constant_ds_rec(expr, fd, cr, flags, t_start, t_end, x, y, min_dt);
                break;
            }
        }
        /*eval_func(t_start, &x, &y, expr, stackpos);
        if ((x == x) && (y == y)) {
            cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
        }*/
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
        uint8_t tflags = 0;
        if (x < window_x0) tflags |= 0x01;
        if (x > window_x1) tflags |= 0x02;
        if (y < window_y0) tflags |= 0x04;
        if (y > window_y1) tflags |= 0x08;

        // [t_start, t_avg]
        draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x0f) | (tflags << 4), t_start, t_avg, x, y, min_dt);
        // [t_avg, t_end]
        draw_function_constant_ds_rec(expr, fd, cr, (flags & 0xf0) | (tflags), t_avg, t_end, x, y, min_dt);
    } else {
        // If neither point is in-bounds, subdivide and iterate
        double t_avg = (t_end + t_start)/2;
        eval_func(t_avg, &x, &y, expr, stackpos);
        if (hypot((x - xi)*xscale, (y - yi)*yscale) < PLOT_MAX_ARC_LENGTH) return;
        //printf("Neither point in bounds on interval %f, %f -> (%f, %f), (%f, %f), %d\n", t_start, t_end, x, y, xi, yi, IN_BOUNDS(x, y));
        uint8_t tflags = 0;
        if (x < window_x0) tflags |= 0x01;
        if (x > window_x1) tflags |= 0x02;
        if (y < window_y0) tflags |= 0x04;
        if (y > window_y1) tflags |= 0x08;
        // [t_start, t_avg]
        draw_function_constant_ds_rec(expr, fd, cr, (flags & 0x0f) | (tflags << 4), t_start, t_avg, x, y, min_dt);
        // [t_avg, t_end]
        draw_function_constant_ds_rec(expr, fd, cr, (flags & 0xf0) | (tflags), t_avg, t_end, x, y, min_dt);
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
        expr->flags |= EXPRESSION_EXPLICIT;
    }
    nfev = 0;
    niev = 0;
    uint8_t flags;
    uint32_t n = (type>>8)/GET_STEP(type);
    //printf("Number of functions to be plotted: %d\n", n);
    expr->cache_size = 0;
    if (expr->special_points) free(expr->special_points);
    expr->special_points = NULL;
    expr->n_special_points = 0;
    int func_pos = fd->n_func;
    for (int i=0; i < fd->n_expr; i++) 
        evaluate_branch(fd->function_list, fd->expression_list[i].func, &func_pos, &stackpos, 1);
    double *old_stack = fd->stack;
    fd->n_stack = stackpos - fd->stack;
    printf("evaluate_branch added %d function blocks\n", func_pos - fd->n_func);
    for (uint32_t i=0; i < n; i++) {
        flags = 0;
        eval_index = i;
        eval_func(t, &xp, &yp, expr, stackpos);
        if (xp < window_x0) flags |= 0x01;
        if (xp > window_x1) flags |= 0x02;
        if (yp < window_y0) flags |= 0x04;
        if (yp > window_y1) flags |= 0x08;
        if (!((xp == xp) && (yp == yp))) flags |= 0x0f;
        eval_func(t_end, &x, &y, expr, stackpos);
        if (x < window_x0) flags |= 0x10;
        if (x > window_x1) flags |= 0x20;
        if (y < window_y0) flags |= 0x40;
        if (y > window_y1) flags |= 0x80;
        if (!((x == x) && (y == y))) flags |= 0xf0;
        //printf("Initial values: (%f, %f) and (%f, %f), flags %02x\n", x, y, xp, yp, flags);
        draw_function_constant_ds_rec(expr, fd, cr, flags, t, t_end, xp, yp, min_dt);
    }
    for (int i=0; i < expr->n_special_points; i++) {
        printf("Special point at %f: (%f, %f)\n", expr->special_points[3*i], expr->special_points[3*i+1], expr->special_points[3*i+2]);
    }
    function *tempfunc;
    for (int i=0; i < fd->n_func; i++) {
        tempfunc = fd->function_list + i;
        if ((tempfunc->oper == func_value) && (tempfunc->first_arg >= fd->function_list + fd->n_func) && (tempfunc->first_arg < fd->function_list + func_pos)) {
            tempfunc = tempfunc->first_arg;
            memcpy(fd->function_list + i, tempfunc, sizeof(function));
        }
    }
    fd->stack = old_stack;
}

void follow_contour(expression *expr, file_data *fd, cairo_t *cr, uint8_t *color, double *area, double e00, double e10, double e11, double e01, uint8_t start_edge) {
    double *stackpos = fd->stack + fd->n_stack;
    double x0 = area[0], x1 = area[1], y0 = area[2], y1 = area[3];

    double bx0 = x0, bx1 = x1, by0 = y0, by1 = y1;
    double be00 = e00, be10 = e10, be11 = e11, be01 = e01;
    double center;
    switch (start_edge) {
        case 0:
            // north edge
            for (uint8_t i=0; i < 3; i++) {
                center = (bx0 + bx1)/2;
                eval_func_2d(&center, &y1, expr, stackpos);
                if (((stackpos[0] >= 0) && (be01 <= 0)) || ((stackpos[0] <= 0) && (be01 >= 0))) {
                    // Contour is between bx0 and center
                    bx1 = center;
                    be11 = stackpos[0];
                } else {
                    // Contour is between center and bx1
                    bx0 = center;
                    be01 = stackpos[0];
                }
                by0 = (by0 + by1)/2;
            }
            eval_func_2d(&bx0, &by0, expr, stackpos);
            be00 = stackpos[0];
            eval_func_2d(&bx1, &by0, expr, stackpos);
            be10 = stackpos[0];
            printf("North edge (%f) converged from [%f, %f] to [%f, %f]\n", y1, x0, x1, bx0, bx1);
            break;
        case 1:
            // east edge
            for (uint8_t i=0; i < 3; i++) {
                center = (by0 + by1)/2;
                eval_func_2d(&x1, &center, expr, stackpos);
                if (((stackpos[0] >= 0) && (be10 <= 0)) || ((stackpos[0] <= 0) && (be10 >= 0))) {
                    // Contour is between by0 and center
                    by1 = center;
                    be11 = stackpos[0];
                } else {
                    // Contour is between center and by1
                    by0 = center;
                    be10 = stackpos[0];
                }
                bx0 = (bx0 + bx1)/2;
            }
            eval_func_2d(&bx0, &by0, expr, stackpos);
            be00 = stackpos[0];
            eval_func_2d(&bx0, &by1, expr, stackpos);
            be01 = stackpos[0];
            printf("East edge (%f) converged from [%f, %f] to [%f, %f]\n", x1, y0, y1, by0, by1);
            break;
        case 2:
            // south edge
            for (uint8_t i=0; i < 3; i++) {
                center = (bx0 + bx1)/2;
                eval_func_2d(&center, &y0, expr, stackpos);
                if (((stackpos[0] >= 0) && (be00 <= 0)) || ((stackpos[0] <= 0) && (be00 >= 0))) {
                    // Contour is between bx0 and center
                    bx1 = center;
                    be10 = stackpos[0];
                } else {
                    // Contour is between center and bx1
                    bx0 = center;
                    be00 = stackpos[0];
                }
                by1 = (by0 + by1)/2;
            }
            eval_func_2d(&bx0, &by1, expr, stackpos);
            be01 = stackpos[0];
            eval_func_2d(&bx1, &by1, expr, stackpos);
            be11 = stackpos[0];
            printf("South edge (%f) converged from [%f, %f] to [%f, %f]\n", y0, x0, x1, bx0, bx1);
            break;
        case 3:
            // west edge
            for (uint8_t i=0; i < 3; i++) {
                center = (by0 + by1)/2;
                eval_func_2d(&x0, &center, expr, stackpos);
                if (((stackpos[0] >= 0) && (be00 <= 0)) || ((stackpos[0] <= 0) && (be00 >= 0))) {
                    // Contour is between by0 and center
                    by1 = center;
                    be01 = stackpos[0];
                } else {
                    // Contour is between center and by1
                    by0 = center;
                    be00 = stackpos[0];
                }
                bx1 = (bx0 + bx1)/2;
            }
            eval_func_2d(&bx1, &by0, expr, stackpos);
            be10 = stackpos[0];
            eval_func_2d(&bx1, &by1, expr, stackpos);
            be11 = stackpos[0];
            printf("West edge (%f) converged from [%f, %f] to [%f, %f]\n", x0, y0, y1, by0, by1);
            break;
        default:
            break;
    }
    double dx = bx1-bx0;
    double dy = by1-by0;
    double npos, epos, spos, wpos;
    uint8_t edges = 0, edge_mask = ~(1<<start_edge);
    for (uint8_t i=0; i < 10; i++) {
        cairo_rectangle(cr, SCALE_XK(bx0), SCALE_YK(by1), xscale*(bx1 - bx0), yscale*(by1 - by0));
        cairo_stroke(cr);
        //return;
        edges = 0;
        if ((be00 != be10) && (((be00 <= 0) && (be10 >= 0)) || ((be00 >= 0) && (be10 <= 0)))) {
            // south edge
            edges |= 0x4;
            spos = bx0 - be00*(bx1 - bx0)/(be10 - be00);
        }
        if ((be10 != be11) && (((be10 <= 0) && (be11 >= 0)) || ((be10 >= 0) && (be11 <= 0)))) {
            // east edge
            edges |= 0x2;
            epos = by0 - be10*(by1 - by0)/(be11 - be10);
        }
        if ((be11 != be01) && (((be11 <= 0) && (be01 >= 0)) || ((be11 >= 0) && (be01 <= 0)))) {
            // north edge
            edges |= 0x1;
            npos = bx0 - be01*(bx0 - bx1)/(be01 - be11);
        }
        if ((be01 != be00) && (((be01 <= 0) && (be00 >= 0)) || ((be01 >= 0) && (be00 <= 0)))) {
            // west edge
            edges |= 0x8;
            wpos = by0 - be00*(by0 - by1)/(be00 - be01);
        }
        double sx0 = SCALE_XK(bx0), sx1 = SCALE_XK(bx1), sy0 = SCALE_YK(by0), sy1 = SCALE_YK(by1);
        double snpos = SCALE_XK(npos), sepos = SCALE_YK(epos), sspos = SCALE_XK(spos), swpos = SCALE_YK(wpos);
        switch (edges) {
            case 0xf:
                // If all four edges are selected, we know that all of the corners have
                // different values from their neighbors. No more than two corners can
                // be zero and they must be on opposite corners. If we sample two adjacent
                // corners, at least one will not be zero
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
                cairo_stroke(cr);
                break;
            case 0xb:
            case 0xe:
            case 0xa:
                cairo_move_to(cr, sx0, swpos);
                cairo_line_to(cr, sx1, sepos);
                cairo_stroke(cr);
                break;
            case 0x3:
                cairo_move_to(cr, sx1, sepos);
                cairo_line_to(cr, snpos, sy1);
                cairo_stroke(cr);
                break;
            case 0x6:
                cairo_move_to(cr, sx1, sepos);
                cairo_line_to(cr, sspos, sy0);
                cairo_stroke(cr);
                break;
            case 0xc:
                cairo_move_to(cr, sx0, swpos);
                cairo_line_to(cr, sspos, sy0);
                cairo_stroke(cr);
                break;
            case 0x9:
                cairo_move_to(cr, sx0, swpos);
                cairo_line_to(cr, snpos, sy1);
                cairo_stroke(cr);
                break;
            default:
                break;
        }
        edges &= edge_mask;
        edge_mask = 0xf;
        printf("[%f, %f], [%f, %f] --> e00 %f, e10 %f, e11 %f, e01 %f --> %x\n", bx0, bx1, by0, by1, be00, be10, be11, be01, edges);
        if (edges & 0x1) {
            // Contour passes through the north edge
            by0 += dy;
            by1 += dy;
            be00 = be01;
            be10 = be11;
            eval_func_2d(&bx0, &by1, expr, stackpos);
            be01 = stackpos[0];
            eval_func_2d(&bx1, &by1, expr, stackpos);
            be11 = stackpos[0];
            edge_mask &= ~0x4;
        } else if (edges & 0x2) {
            // Contour passes through the east edge
            bx0 += dx;
            bx1 += dx;
            be00 = be10;
            be01 = be11;
            eval_func_2d(&bx1, &by0, expr, stackpos);
            be10 = stackpos[0];
            eval_func_2d(&bx1, &by1, expr, stackpos);
            be11 = stackpos[0];
            edge_mask &= ~0x8;
        } else if (edges & 0x4) {
            // Contour passes through the south edge
            by0 -= dy;
            by1 -= dy;
            be01 = be00;
            be11 = be10;
            eval_func_2d(&bx0, &by0, expr, stackpos);
            be00 = stackpos[0];
            eval_func_2d(&bx1, &by0, expr, stackpos);
            be10 = stackpos[0];
            edge_mask &= ~0x1;
        } else if (edges & 0x8) {
            // Contour passes through the west edge
            bx0 -= dx;
            bx1 -= dx;
            be10 = be00;
            be11 = be01;
            eval_func_2d(&bx0, &by0, expr, stackpos);
            be00 = stackpos[0];
            eval_func_2d(&bx0, &by1, expr, stackpos);
            be01 = stackpos[0];
            edge_mask &= ~0x2;
        }
    }
    cairo_rectangle(cr, SCALE_XK(bx0), SCALE_YK(by1), xscale*(bx1 - bx0), yscale*(by1 - by0));
    cairo_stroke(cr);
}

void draw_implicit_rec(expression *expr, file_data *fd, cairo_t *cr, uint8_t *color, double *area, int divisions, double *vbuf, double *hbuf, uint8_t buf_valid) {
    /*fd->variable_list[0].pointer = area;
    fd->variable_list[0].type = 1<<8;
    fd->variable_list[1].pointer = area+2;
    fd->variable_list[1].type = 1<<8;*/
    double *stackpos = fd->stack + fd->n_stack;
    double *lstackpos = fd->lstack;
    eval_inter_2d(area, area+2, expr, stackpos, lstackpos);
    double z0 = lstackpos[0], z1 = stackpos[0];
    uint32_t bufsize = (1<<(PLOT_IMPLICIT_HARDMAX - divisions))+1;
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
            enclosed_area += (area[1] - area[0])*(area[3] - area[2]);
        }
#endif
    }
    if ((stackpos[0] < 0) || (lstackpos[0] > 0)) {
        // No contours can exist in the given area
        //printf("    No contours found in ([%f, %f], [%f, %f])\n", area[0], area[1], area[2], area[3]);
        //cairo_rectangle(cr, SCALE_XK(area[0]), SCALE_YK(area[3]), xscale*(area[1] - area[0]), yscale*(area[3] - area[2]));
        //cairo_stroke(cr);
        hbuf[0] = vbuf[bufsize-1];
        vbuf[0] = hbuf[bufsize-1];
        for (uint32_t i=1; i<bufsize; i++) {
            hbuf[i] = NAN;
            vbuf[i] = NAN;
        }
        return;
    }
    if ((divisions == PLOT_IMPLICIT_MAXDEPTH) || (divisions >= PLOT_IMPLICIT_HARDMAX)) {
        // Maximum depth reached
        //printf("    maximum depth reached on ([%f, %f], [%f, %f]) --> [%f, %f]\n", area[0], area[1], area[2], area[3], lstackpos[0], stackpos[0]);
        //cairo_rectangle(cr, SCALE_XK(area[0]), SCALE_YK(area[3]), xscale*(area[1] - area[0]), yscale*(area[3] - area[2]));
        //cairo_stroke(cr);
        function *arg1 = expr->func->first_arg;
        function *arg2 = arg1->next_arg;
        double x0 = area[0], x1 = area[1], y0 = area[2], y1 = area[3];
        double e00, e10, e11, e01;
        if ((buf_valid & 0x01) && (vbuf[0] == vbuf[0])) e00 = vbuf[0];
        else if ((buf_valid & 0x02) && (hbuf[0] == hbuf[0])) e00 = hbuf[0];
        else {
            eval_func_2d(&x0, &y0, expr, stackpos); 
            e00 = stackpos[0];
        }
        if (hbuf[bufsize-1] == hbuf[bufsize-1]) e10 = hbuf[bufsize-1];
        else {
            eval_func_2d(&x1, &y0, expr, stackpos); 
            e10 = stackpos[0];
        }
        if (vbuf[bufsize-1] == vbuf[bufsize-1]) {
            e01 = vbuf[bufsize-1];
        } else {
            eval_func_2d(&x0, &y1, expr, stackpos); 
            e01 = stackpos[0];
        }
        eval_func_2d(&x1, &y1, expr, stackpos); 
        e11 = stackpos[0];

        if ((e00 > 0) && (e10 > 0) && (e11 > 0) && (e01 > 0)) {
            // Interval calculation was wrong
#ifdef PLOT_USE_INEQUALITY
            if (ineq) {
                SET_COLOR_OPACITY(cr, color, 0.5);
                cairo_rectangle(cr, SCALE_XK(area[0]), SCALE_YK(area[3]), xscale*(area[1] - area[0]), yscale*(area[3] - area[2]));
                cairo_fill(cr);
                SET_COLOR_OPACITY(cr, color, 1);
                enclosed_area += (area[1] - area[0])*(area[3] - area[2]);
            }
#endif
            vbuf[0] = e10;
            vbuf[bufsize-1] = e11;
            hbuf[0] = e01;
            hbuf[bufsize-1] = e11;
            for (int i=1; i < bufsize - 1; i++) {
                hbuf[i] = NAN;
                vbuf[i] = NAN;
            }
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
        /*if (edges & 0x1) follow_contour(expr, fd, cr, color, area, e00, e10, e11, e01, 0);
        else if (edges & 0x2) follow_contour(expr, fd, cr, color, area, e00, e10, e11, e01, 1);
        else if (edges & 0x4) follow_contour(expr, fd, cr, color, area, e00, e10, e11, e01, 2);
        else if (edges & 0x8) follow_contour(expr, fd, cr, color, area, e00, e10, e11, e01, 3);*/
        //if (divisions == PLOT_IMPLICIT_MAXDEPTH) printf("([%f, %f], [%f, %f]) --> edgex %x, e00 %f, e10 %f, e11 %f, e01 %f\n", area[0], area[1], area[2], area[3], edges, e00, e10, e11, e01);
        if (((singular) || (edges)) && (divisions < PLOT_IMPLICIT_HARDMAX)) {
            //printf("    Unresolved detail in ([%f, %f], [%f, %f])\n", area[0], area[1], area[2], area[3]);
            // Subdivide
            double x0 = area[0], x1 = area[1], y0 = area[2], y1 = area[3];
            double xm = (x0 + x1)/2, ym = (y0 + y1)/2;
            //printf("subdividing ([%f, %f], [%f, %f]), %f, %f\n", x0, x1, y0, y1, xm, ym);
            // Bottom left box
            uint8_t temp_buf_valid;
            area[1] = xm; area[3] = ym;
            draw_implicit_rec(expr, fd, cr, color, area, divisions+1, vbuf, hbuf, buf_valid);
            // Bottom right box
            area[0] = xm; area[1] = x1; area[2] = y0; area[3] = ym;
            draw_implicit_rec(expr, fd, cr, color, area, divisions+1, vbuf, hbuf+(bufsize>>1), 0x01);
            // Top left box
            area[0] = x0; area[1] = xm; area[2] = ym; area[3] = y1;
            draw_implicit_rec(expr, fd, cr, color, area, divisions+1, vbuf+(bufsize>>1), hbuf, 0x02);
            // Top right box
            area[0] = xm; area[1] = x1; area[2] = ym; area[3] = y1;
            draw_implicit_rec(expr, fd, cr, color, area, divisions+1, vbuf+(bufsize>>1), hbuf+(bufsize>>1), 0x01);
            return;
        } else {
            if ((!singular) && (!edges)) {
                printf("Skipping box ([%f, %f], [%f, %f]) with corners %e, %e, %e, %e in (%e, %e)\n", x0, x1, y0, y1, e00, e10, e11, e01, z0, z1);
            }
            vbuf[0] = e10;
            vbuf[bufsize-1] = e11;
            hbuf[0] = e01;
            hbuf[bufsize-1] = e11;
            for (int i=1; i < bufsize - 1; i++) {
                hbuf[i] = NAN;
                vbuf[i] = NAN;
            }
        }
        // If two edges are selected, then draw the line between those two edges. If three 
        // edges are selected, then one of the corners must be zero so the contour must go
        // through one of the corners. If four edges are selected, there are two contour
        // lines, each going through opposite sides of the area.
        double sx0 = SCALE_XK(x0), sx1 = SCALE_XK(x1), sy0 = SCALE_YK(y0), sy1 = SCALE_YK(y1);
        double snpos = SCALE_XK(npos), sepos = SCALE_YK(epos), sspos = SCALE_XK(spos), swpos = SCALE_YK(wpos);
        // If any of the corners is undefined, draw nothing
        if ((e00 != e00) || (e01 != e01) || (e10 != e10) || (e11 != e11)) return;
        if (isinf(e00) || isinf(e01) || isinf(e10) || isinf(e11)) return;
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
                        cairo_line_to(cr, sx0, swpos);
                        cairo_line_to(cr, sx1, sepos);
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
                        enclosed_area += ((spos + npos)/2 - x0) * (y1 - y0);
                    } else {
                        cairo_line_to(cr, sx1, sy1);
                        cairo_line_to(cr, sx1, sy0);
                        enclosed_area += (x1 - (spos + npos)/2) * (y1 - y0);
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
                        enclosed_area += ((epos + wpos)/2 - y0) * (x1 - x0);
                    } else {
                        cairo_line_to(cr, sx1, sy1);
                        cairo_line_to(cr, sx0, sy1);
                        enclosed_area += (y1 - (epos + wpos)/2) * (x1 - x0);
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
                        enclosed_area += (x1 - npos)*(y1 - epos)/2;
                    } else {
                        cairo_line_to(cr, sx0, sy1);
                        cairo_line_to(cr, sx0, sy0);
                        cairo_line_to(cr, sx1, sy0);
                        enclosed_area += (y1 - y0)*(x1 - x0) - (x1 - npos)*(y1 - epos)/2;
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
                        enclosed_area += (spos - x1)*(epos - y0)/2;
                    } else {
                        cairo_line_to(cr, sx0, sy0);
                        cairo_line_to(cr, sx0, sy1);
                        cairo_line_to(cr, sx1, sy1);
                        enclosed_area += (y1 - y0)*(x1 - x0) - (spos - x1)*(epos - y0)/2;
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
                        enclosed_area += (spos - x0)*(wpos - y0)/2;
                    } else {
                        cairo_line_to(cr, sx1, sy0);
                        cairo_line_to(cr, sx1, sy1);
                        cairo_line_to(cr, sx0, sy1);
                        enclosed_area += (y1 - y0)*(x1 - x0) - (spos - x0)*(wpos - y0)/2;
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
                        enclosed_area += (npos - x0)*(y1 - wpos)/2;
                    } else {
                        cairo_line_to(cr, sx1, sy1);
                        cairo_line_to(cr, sx1, sy0);
                        cairo_line_to(cr, sx0, sy0);
                        enclosed_area += (y1 - y0)*(x1 - x0) - (npos - x0)*(y1 - wpos)/2;
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
        // Bottom left box
        uint8_t temp_buf_valid;
        area[1] = xm; area[3] = ym;
        draw_implicit_rec(expr, fd, cr, color, area, divisions+1, vbuf, hbuf, buf_valid);
        // Bottom right box
        area[0] = xm; area[1] = x1; area[2] = y0; area[3] = ym;
        draw_implicit_rec(expr, fd, cr, color, area, divisions+1, vbuf, hbuf+(bufsize>>1), 0x01);
        // Top left box
        area[0] = x0; area[1] = xm; area[2] = ym; area[3] = y1;
        draw_implicit_rec(expr, fd, cr, color, area, divisions+1, vbuf+(bufsize>>1), hbuf, 0x02);
        // Top right box
        area[0] = xm; area[1] = x1; area[2] = ym; area[3] = y1;
        draw_implicit_rec(expr, fd, cr, color, area, divisions+1, vbuf+(bufsize>>1), hbuf+(bufsize>>1), 0x01);
    }
}

void draw_implicit(expression *expr, file_data *fd, uint8_t *color, cairo_t *cr) {
    if ((expr->color_pointer) && (expr->color_length >= 3)) {
        color[0] = (uint8_t)(expr->color_value[0]);
        color[1] = (uint8_t)(expr->color_value[1]);
        color[2] = (uint8_t)(expr->color_value[2]);
    }
    SET_COLOR(cr, color);
    double temp[4] = {window_x0, window_x1, window_y0, window_y1};
    uint32_t (*old_oper)(void*, double*) = expr->func->oper;
    if (expr->func->oper == func_equals) expr->func->oper = func_sub;
    else if (expr->func->oper == func_compare_single) expr->func->oper = func_compare_sub_single;
    else if (expr->func->oper == func_compare) expr->func->oper = func_compare_sub;
    double *vbuf = malloc(sizeof(double)*((1<<PLOT_IMPLICIT_HARDMAX)+1));
    double *hbuf = malloc(sizeof(double)*((1<<PLOT_IMPLICIT_HARDMAX)+1));
    for (uint32_t i=0; i <= (1<<PLOT_IMPLICIT_HARDMAX); i++) {
        vbuf[i] = NAN;
        hbuf[i] = NAN;
    }
    enclosed_area = 0;
    draw_implicit_rec(expr, fd, cr, color, temp, 0, vbuf, hbuf, 0x00);
    printf("Enclosed area is %f\n", enclosed_area);
    free(vbuf);
    free(hbuf);
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
    g_print ("Clicked at %f, %f, state %02x\n", event->x, event->y, event->state);
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
    click_state |= 0x01;
    if ((event->state & GDK_SHIFT_MASK) && (fabs(event->x - SCALE_XK(0)) < AXIS_MARGIN)) {
        click_state |= 0x02;
    } else if ((event->state & GDK_SHIFT_MASK) && (fabs(event->y - SCALE_YK(0)) < AXIS_MARGIN)) {
        click_state |= 0x04;
    } else if (event->state & GDK_SHIFT_MASK) click_state |= 0x06;
    return TRUE;
}

static gboolean button_release_callback (GtkWidget *event_box, GdkEventButton *event, gpointer data_pointer) {
    file_data *fd = (file_data*)(data_pointer);
    printf("Clearing overlay\n");
    reset_overlay(fd);
    click_x = event->x;
    click_y = event->y;
    click_state &= 0xf8;
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
    cairo_text_extents_t extents;
    char temp[50];
    int16_t log, decade;
    int8_t base;
    double major_xticksize, minor_xticksize, major_yticksize, minor_yticksize;
    clock_t t3, t4;

    log = round(3*log10(TICK_SIZE / xscale));
    base = log%3;
    base = (base<0 ? 3+base : base);
    decade = (log - base)/3;
    int16_t digitsx = (decade < 0 ? -decade : 0);
    base++;
    if (base == 3) base = 5;
    major_xticksize = base;
    while (decade > 0) {decade--; major_xticksize *= 10;}
    while (decade < 0) {decade++; major_xticksize /= 10;}
    if (base == 2) minor_xticksize = major_xticksize / 2;
    else minor_xticksize = major_xticksize / 5;

    log = round(3*log10(TICK_SIZE / yscale));
    base = log%3;
    base = (base<0 ? 3+base : base);
    decade = (log - base)/3;
    int16_t digitsy = (decade < 0 ? -decade : 0);
    base++;
    if (base == 3) base = 5;
    major_yticksize = base;
    while (decade > 0) {decade--; major_yticksize *= 10;}
    while (decade < 0) {decade++; major_yticksize /= 10;}
    if (base == 2) minor_yticksize = major_yticksize / 2;
    else minor_yticksize = major_yticksize / 5;
    
    cairo_set_source_rgb(cr, COLOR_MINOR_GRID/256.0, COLOR_MINOR_GRID/256.0, COLOR_MINOR_GRID/256.0);
    double x0_scaled = minor_xticksize*floor(window_x0/minor_xticksize);
    for (double tick=x0_scaled; tick <= window_x1; tick+=minor_xticksize) {
        cairo_move_to(cr, SCALE_XK(tick), 0);
        cairo_line_to(cr, SCALE_XK(tick), HEIGHT);
    }
    double y0_scaled = minor_yticksize*floor(window_y0/minor_yticksize);
    for (double tick=y0_scaled; tick <= window_y1; tick+=minor_yticksize) {
        cairo_move_to(cr, 0, SCALE_YK(tick));
        cairo_line_to(cr, WIDTH, SCALE_YK(tick));
    }
    cairo_stroke(cr);
    
    cairo_set_source_rgb(cr, COLOR_MAJOR_GRID/256.0, COLOR_MAJOR_GRID/256.0, COLOR_MAJOR_GRID/256.0);
    x0_scaled = major_xticksize*floor(window_x0/major_xticksize);
    for (double tick=x0_scaled; tick <= window_x1; tick+=major_xticksize) {
        cairo_move_to(cr, SCALE_XK(tick), 0);
        cairo_line_to(cr, SCALE_XK(tick), HEIGHT);

        cairo_select_font_face(cr, "Verdana", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 12);
        sprintf(temp, "%.*f", digitsx, tick);
        cairo_text_extents(cr, temp, &extents);
        double ypos = SCALE_YK(0)+extents.height+AXIS_DISTANCE;
        if (ypos >= HEIGHT) ypos = HEIGHT - AXIS_DISTANCE;
        else if (window_y1 <= 0) ypos = extents.height + AXIS_DISTANCE;
        cairo_move_to(cr, SCALE_XK(tick)-extents.width/2, ypos);
        cairo_show_text(cr, temp);
    }
    y0_scaled = major_yticksize*floor(window_y0/major_yticksize);
    for (double tick=y0_scaled; tick <= window_y1; tick+=major_yticksize) {
        cairo_move_to(cr, 0, SCALE_YK(tick));
        cairo_line_to(cr, WIDTH, SCALE_YK(tick));

        cairo_select_font_face(cr, "Verdana", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 12);
        sprintf(temp, "%.*f", digitsy, tick);
        cairo_text_extents(cr, temp, &extents);
        double xpos = SCALE_XK(0)-extents.width-AXIS_DISTANCE;
        if (xpos < 0) xpos = AXIS_DISTANCE;
        else if (window_x1 <= 0) xpos = WIDTH-extents.width-AXIS_DISTANCE;
        cairo_move_to(cr, xpos, SCALE_YK(tick)+extents.height/2);
        cairo_show_text(cr, temp);
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
            expr->color_value = stack+n_stack;
            expr->color_length = color_len;
            fd->n_stack += color_len;
            n_stack += color_len;
        } else {
            color_len = 0;
            expr->color_length = 0;
            expr->color_value = NULL;
        }
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
            expr->flags |= EXPRESSION_EVALUATE;
            if (!from) from = expr;
            if ((((expr->var->type)>>8) != 0) && (expr->var->pointer)) free(expr->var->pointer);
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
    printf("drawingarea size is (%d, %d)\n", gtk_widget_get_allocated_width(widget), gtk_widget_get_allocated_height(widget));
    da_width = gtk_widget_get_allocated_width(widget);
    da_height = gtk_widget_get_allocated_height(widget);
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
    if ((click_state == 1) && (fd->use_overlay) && (fd->click_expr)) {
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
        // Draw special points
        SET_COLOR_HEX(cr, COLOR_SPECIAL_POINT);
        for (uint32_t i=0; i < fd->click_expr->n_special_points; i++) {
            cairo_arc(cr, SCALE_XK(fd->click_expr->special_points[3*i+1]), SCALE_YK(fd->click_expr->special_points[3*i+2]), POINT_SIZE, 0, 2*G_PI);
            cairo_fill(cr);
        }
        // Draw rectangle
        if (!status) draw_labeled_point(cr, xf, yf);
        cairo_destroy(cr);
        gtk_widget_queue_draw(fd->drawing_area);
    } else if (click_state & 0x01) {
        double center_x = 0, center_y = 0;
        if (window_x0 > 0) center_x = window_x0;
        if (window_x1 < 0) center_x = window_x1;
        if (window_y0 > 0) center_y = window_y0;
        if (window_y1 < 0) center_y = window_y1;
        if (event->state & GDK_SHIFT_MASK) {
            // Rescale uniformly
            // Scale factor is r_new / r_old
            if (!(click_state & 0x02)) click_y = event->y;
            if (!(click_state & 0x04)) click_x = event->x;
            double scale = hypot(SCALE_XK(center_x) - event->x, SCALE_YK(center_y) - event->y)/hypot(SCALE_XK(center_x) - click_x, SCALE_YK(center_y) - click_y);
            if (click_state & 0x04) {
                window_x0 = center_x + (window_x0 - center_x)/scale;
                window_x1 = center_x + (window_x1 - center_x)/scale;
            }
            if (click_state & 0x02) {
                window_y0 = center_y + (window_y0 - center_y)/scale;
                window_y1 = center_y + (window_y1 - center_y)/scale;
            }
        } else {
            // Shift key is not pressed
            double dx = (event->x - click_x)/xscale;
            double dy = (click_y - event->y)/yscale;
            window_x0 -= dx; window_x1 -= dx; window_y0 -= dy; window_y1 -= dy;
            uint8_t color[4] = {255, 255, 255, 0};
        }
        click_x = event->x;
        click_y = event->y;
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
    printf("Event is %d, state %08x\n", event->keyval, event->state);
    if (event->keyval == 's') {
        if (event->state & GDK_CONTROL_MASK) {
            cairo_surface_write_to_png(fd->surface, "/tmp/plotter_export.png");
        } else {
            run_ticker = 1;
            if (ticker_target >= 0) g_timeout_add(ticker_step, timeout_callback, data_pointer);
        }
    } else if (event->keyval == 'e') {
        run_ticker = 0;
    } else if (event->keyval == 'i') {
        treeview_activate(fd);
        printf("inspecting\n");
    } else if (event->keyval == GDK_KEY_Home) {
        double w = da_width / INIT_SCALE, h = da_height / INIT_SCALE;
        window_x0 = -w/2;
        window_x1 = w/2;
        window_y0 = -h/2;
        window_y1 = h/2;
        xscale = INIT_SCALE;
        yscale = INIT_SCALE;
        redraw_all(fd);
        gtk_widget_queue_draw(fd->drawing_area);
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

void terminate(int sig) {
    // Deallocate any memory
    printf("Closing\n");
    variable *variable_list = fd.variable_list;
    function *function_list = fd.function_list;
    // Free memory allocated to variables
    for (int i=4; i < fd.n_var; i++) {
        if ((variable_list[i].pointer) && ((variable_list[i].type>>8) != 0) && !(variable_list[i].flags & (VARIABLE_ARGUMENT | VARIABLE_FUNCTION | VARIABLE_ACTION))) {
            free(variable_list[i].pointer);
        }
    }
    // Free memory allocated for integration data
    for (int i=0; i < fd.n_func; i++) {
        if ((function_list[i].oper == func_integrate_gsl) && (function_list[i].value)) free(function_list[i].value);
    }
    exit(0);
}

int main (int argc, char **argv) {
    signal(SIGINT, terminate);

    memset(variable_list, 0, 10*sizeof(variable));
    memset(expression_list, 0, 10*sizeof(expression));
    printf("First stack positions %p, %p\n", stack, stack+1);
    variable_list[0] = new_variable("x", 0, VARIABLE_IN_SCOPE | VARIABLE_INTERVAL | VARIABLE_NOT_FIXED, NULL);
    variable_list[1] = new_variable("y", 0, VARIABLE_IN_SCOPE | VARIABLE_INTERVAL | VARIABLE_NOT_FIXED, NULL);
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


    app = gtk_application_new("org.gtk.example", G_APPLICATION_FLAGS_NONE);
    g_signal_connect(app, "activate", G_CALLBACK (activate), &fd);
    for (int i=4; i < argc; i++) argv[i-3] = argv[i];
    status = g_application_run(G_APPLICATION (app), argc-3, argv);
    g_object_unref (app);

    terminate(0);
}
