#include <gtk/gtk.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "parse.h"
#include "functions.h"
#include "linalg_functions.h"
#include "config.h"

#define SCALE_X(v) ((uint16_t)(250*(1+v)))
#define SCALE_Y(v) ((uint16_t)(250*(1-v)))

#define SCALE_XF(v) ((int64_t)(250L*(1<<16)*(1+v)))
#define SCALE_YF(v) ((int64_t)(250L*(1<<16)*(1-v)))

#define SCALE_XK(v) ((v - window_x0)*xscale)
#define SCALE_YK(v) ((window_y1 - v)*yscale)

#define CLIP_WIDTH(v) (v<0 ? 0 : (v>WIDTH ? WIDTH : v))
#define CLIP_HEIGHT(v) (v<0 ? 0 : (v>HEIGHT ? HEIGHT : v))
#define IN_BOUNDS(x, y) ((window_x0 <= x) && (x < window_x1) && (window_y0 <= y) && (y < window_y1))
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

function function_list[500];
variable variable_list[50];
expression expression_list[100];
double stack[65536];
char stringbuf[500];

uint32_t n_expr;
expression *top_expr;

clock_t t1, t2;

uint8_t click_state = 0;
double click_x;
double click_y;
int ticker_target, ticker_step;
uint8_t run_ticker;

cairo_matrix_t transform_matrix;


void eval_func(double t, double *x, double *y, function *func, double *stackpos) {
    variable_list[0].pointer = &t;
    variable_list[0].type = 1<<8;
    func->oper(func, stackpos);
    *y = stackpos[0];
    *x = t;
    nfev++;
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
    double x, y, xp, yp;
    eval_func(window_x0, &xp, &yp, func, stackpos);
    cairo_move_to(cr, SCALE_XK(xp), SCALE_YK(yp));
    double minstep_x = (window_x1 - window_x0)/WIDTH;
    for (double t = window_x0+minstep_x; t < window_x1; t+=minstep_x) {
        eval_func(t, &x, &y, func, stackpos);
        if ((x+1 > x) && (y+1 > y) && (xp+1 > xp) && (yp+1 > yp) && (IN_BOUNDS(x, y) || IN_BOUNDS(xp, yp))) cairo_line_to(cr, SCALE_XK(x), SCALE_YK(y));
        else if ((x+1 > x) && (y+1 > y)) cairo_move_to(cr, SCALE_XK(x), SCALE_YK(y));
        xp = x; yp = y;
    }
    cairo_stroke(cr);
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
    for (int i=0; i < n_expr; i++) {
        if (expression_list[i].flags & EXPRESSION_PLOTTABLE) {
            if ((expression_list[i].flags & EXPRESSION_FIXED) && ((expression_list[i].value_type & TYPE_MASK) == TYPE_POINT)) {
                int len = (expression_list[i].value_type) >> 8;
                double *ptr = expression_list[i].value;
                SET_COLOR(cr, expression_list[i].color);
                for (int p=0; p < len; p+=2) {
                    pt_x = SCALE_XK(ptr[p]);
                    pt_y = SCALE_YK(ptr[p+1]);
                    cairo_arc(cr, pt_x, pt_y, POINT_SIZE, 0, 2*G_PI);
                    cairo_fill(cr);
                }
            } else if ((expression_list[i].flags & EXPRESSION_FIXED) && ((expression_list[i].value_type & TYPE_MASK) == TYPE_POLYGON)) {
                int len = (expression_list[i].value_type) >> 8;
                double *ptr = expression_list[i].value;
                double x_temp, y_temp;
                if (i != 19) SET_COLOR(cr, expression_list[i].color);
                uint32_t color_pos = 0;
                for (int p=0; p < len; p += 2*MAX_POLYGON_SIZE) {
                    pt_x = SCALE_XK(ptr[p]);
                    pt_y = SCALE_YK(ptr[p+1]);
                    if (i == 19) SET_COLOR(cr, (expression_list[20].value+color_pos));
                    color_pos += 3;
                    cairo_move_to(cr, pt_x, pt_y);
                    for (int k=1; k < MAX_POLYGON_SIZE; k++) {
                        x_temp = ptr[p+2*k]; y_temp = ptr[p+2*k+1];
                        if ((x_temp+1 > x_temp) && (y_temp+1 > y_temp)) {
                            pt_x1 = SCALE_XK(x_temp);
                            pt_y1 = SCALE_YK(y_temp);
                            //cairo_move_to(cr, pt_x, pt_y);
                            cairo_line_to(cr, pt_x1, pt_y1);
                            pt_x = pt_x1;
                            pt_y = pt_y1;
                        }
                    }
                    //cairo_move_to(cr, pt_x, pt_y);
                    //cairo_line_to(cr, SCALE_XK(ptr[p]), SCALE_YK(ptr[p+1]));
                    cairo_close_path(cr);
                    cairo_stroke_preserve(cr);
                    cairo_fill(cr);
                }
                //cairo_stroke(cr);
            } else {
                draw_function_slow(expression_list[i].func, stack+60, expression_list[i].color, cr);
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
        gtk_widget_queue_draw(data_pointer);
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
        gtk_widget_queue_draw(data_pointer);
    }
    return TRUE;
}

static gboolean timeout_callback(gpointer data_pointer) {
    // variable 5 is alpha, expression 3 is alpha definition
    printf("timeout_callback\n");
    clock_t t3 = clock();
    expression_list[ticker_target].func->oper(expression_list[ticker_target].func, stack+60);
    expression *expr = top_expr;
    expression *from = NULL;
    while (expr) {
        if ((expr->var) && (expr->var->new_pointer)) {
            printf("expression %p (offset %ld) has changed, expr->var %p\n", expr, expr - expression_list, expr->var);
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
        evaluate_from(expression_list, n_expr, from, stack+60);
        gtk_widget_queue_draw(data_pointer);
        clock_t t4 = clock();
        printf("Total time: %luus\n", t4-t3);
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

    window = gtk_application_window_new (app);
    event_box = gtk_event_box_new();
    gtk_container_add(GTK_CONTAINER(window), event_box);

    gtk_widget_add_events(event_box, GDK_BUTTON_PRESS_MASK);
    gtk_widget_add_events(event_box, GDK_SCROLL_MASK);
    gtk_widget_add_events(event_box, GDK_POINTER_MOTION_MASK);
    gtk_widget_add_events(event_box, GDK_BUTTON_RELEASE_MASK);

    GtkWidget *drawing_area = gtk_drawing_area_new();
    gtk_widget_set_size_request(drawing_area, WIDTH, HEIGHT);
    g_signal_connect(G_OBJECT(event_box), "button_press_event", G_CALLBACK(button_press_callback), drawing_area);
    g_signal_connect(G_OBJECT(event_box), "button_release_event", G_CALLBACK(button_release_callback), drawing_area);
    g_signal_connect(G_OBJECT(event_box), "scroll_event", G_CALLBACK(scroll_callback), drawing_area);
    g_signal_connect(G_OBJECT(event_box), "motion-notify-event", G_CALLBACK(motion_callback), drawing_area);
    g_signal_connect(G_OBJECT(drawing_area), "draw", G_CALLBACK(redraw_all), drawing_area);
    g_signal_connect(window, "key-press-event", G_CALLBACK(keypress_callback), drawing_area);
    gtk_container_add(GTK_CONTAINER(event_box), drawing_area);

    gtk_widget_show_all(window);
}

int main (int argc, char **argv) {
    memset(variable_list, 0, 10*sizeof(variable));
    memset(expression_list, 0, 10*sizeof(expression));
    printf("First stack positions %p, %p\n", stack, stack+1);
    variable_list[0] = new_variable("x", 0, VARIABLE_IN_SCOPE, NULL);
    variable_list[1] = new_variable("y", 0, VARIABLE_IN_SCOPE, NULL);
    //parse_latex("\\frac{\\arctan\\left(\\frac{0.1}{x}\\right)}{2}", function_list);
    //int stack_size = 0;
    //parse_latex("2x+x^{4-x}\\cos\\left(4+x\\right)-\\left(7+3\\right)", function_list, stack, variable_list, &stack_size);
    //parse_latex("2\\left(\\sin\\left(5x\\right)-\\sin\\left(3x\\right)\\right)", function_list, stack, variable_list);
    //parse_latex("\\cos\\left(16x\\right)", function_list, stack, variable_list);
    //parse_latex("N_{c}\\left[A_{pindex}\\left(N_{f},G\\left(N_{s0}\\left[A_{csnnn}\\right],N_{s1}\\left[A_{csnnn}\\right],\\left[1...\\operatorname{length}\\left(A_{csnnn}\\right)\\right]-A_{csnn}\\left[A_{csnnn}\\right]-1,N_{st}\\left[A_{csnnn}\\right]\\right)\\right)\\right]", function_list, stack);
    
    uint32_t n_func = 0;
    uint32_t n_var = 0;
    n_expr = load_file(argv[1], expression_list);
    top_expr = parse_file(function_list, stack, variable_list, stringbuf, expression_list, &n_func, &n_var, n_expr);
    printf("Parsing completed. %d function blocks, %d variables, %d expressions\n", n_func, n_var, n_expr);
    ticker_target = (int)parse_double(argv[2]);
    ticker_step = (int)parse_double(argv[3]);
    printf("Expression %d will be evaluated every %d milliseconds\n", ticker_target, ticker_step);

    uint32_t type;

    GtkApplication *app;
    int status;


    app = gtk_application_new("org.gtk.example", G_APPLICATION_FLAGS_NONE);
    g_signal_connect(app, "activate", G_CALLBACK (activate), NULL);
    for (int i=4; i < argc; i++) argv[i-3] = argv[i];
    status = g_application_run(G_APPLICATION (app), argc-3, argv);
    g_object_unref (app);

    return status;
}
