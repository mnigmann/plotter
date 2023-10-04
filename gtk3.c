#include <gtk/gtk.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "parse.h"
#include "functions.h"
#include "lines.h"

#define WIDTH 500
#define HEIGHT 500

#define DT_DERIV 1e-4
#define RADIUS_SCALE 1
#define MAX_ARC_LENGTH 0.1

#define PI 3.14159265359
#define SCROLL_AMOUNT 1.1

#define SCALE_X(v) ((uint16_t)(250*(1+v)))
#define SCALE_Y(v) ((uint16_t)(250*(1-v)))

#define SCALE_XF(v) ((int64_t)(250L*(1<<16)*(1+v)))
#define SCALE_YF(v) ((int64_t)(250L*(1<<16)*(1-v)))

#define SCALE_XK(v) ((int64_t)((v - window_x0)*xscale))
#define SCALE_YK(v) ((int64_t)((window_y1 - v)*yscale))

#define CLIP_WIDTH(v) (v<0 ? 0 : (v>WIDTH ? WIDTH : v))
#define CLIP_HEIGHT(v) (v<0 ? 0 : (v>HEIGHT ? HEIGHT : v))

#define TOIDX(x, y) (3*(x + y*WIDTH))

guchar data[3*WIDTH*HEIGHT];
GdkPixbuf *pixbuf;
double window_x0 = 0.0;
double window_y0 = 0.0;
double window_x1 = 10;
double window_y1 = 10;
double xscale;
double yscale;
uint16_t nfev = 0;

function function_list[100];
variable variable_list[30];
expression expression_list[100];
double stack[2000];
char stringbuf[500];

uint32_t n_expr;

clock_t t1, t2;

uint8_t click_state = 0;
double click_x;
double click_y;


void draw_line(uint16_t x0, uint16_t y0, uint16_t x1, uint16_t y1) {
    if ((x0 < 0) || (x0 >= WIDTH) || (y0 < 0) || (y0 >= HEIGHT) || (x1 < 0) || (x1 >= WIDTH) || (y1 < 0) || (y1 >= HEIGHT)) return;
    //data[TOIDX(x0, y0)+1] = 255;
    //data[TOIDX(x1, y1)+1] = 255;
    uint16_t xmin, xmax, ymin, ymax, ystart;
    if (x1 > x0) {
        xmin = x0;
        ystart = y0;
        xmax = x1;
    } else {
        xmin = x1;
        ystart = y1;
        xmax = x0;
    }
    if (y1 > y0) {
        ymin = y0;
        ymax = y1;
    } else {
        ymin = y1;
        ymax = y0;
    }
    uint16_t a = ymax - ymin;
    uint16_t b = xmax - xmin;
    uint32_t pos = 3*(xmin + ystart*WIDTH);
    if ((x1 > x0) != (y1 > y0)) {
        // Ascending line
        if (b >= a) {
            // Line is more horizontal
            int16_t f = -b;
            for (uint16_t i=xmin; i <= xmax; i++) {
                data[pos+2] = 255;
                pos+=3;
                f += 2*a;
                if (f > 0) {
                    pos -= 3*WIDTH;
                    f -= 2*b;
                }
            }
        } else {
            // Line is more vertical
            int16_t f = -a;
            for (uint16_t i=ymin; i <= ymax; i++) {
                data[pos+2] = 255;
                pos-=3*WIDTH;
                f += 2*b;
                if (f > 0) {
                    pos += 3;
                    f -= 2*a;
                }
            }
        }
    } else {
        // Descending line
        if (b >= a) {
            // Line is more horizontal
            int16_t f = -b;
            for (uint16_t i=xmin; i <= xmax; i++) {
                data[pos+2] = 255;
                pos+=3;
                f += 2*a;
                if (f > 0) {
                    pos += 3*WIDTH;
                    f -= 2*b;
                }
            }
        } else {
            // Line is more vertical
            int16_t f = -a;
            for (uint16_t i=ymin; i <= ymax; i++) {
                data[pos+2] = 255;
                pos+=3*WIDTH;
                f += 2*b;
                if (f > 0) {
                    pos += 3;
                    f -= 2*a;
                }
            }
        }
    }
}

void eval_func(double t, double *x, double *y, function *func, double *stackpos) {
    variable_list[0].pointer = &t;
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

void draw_function(function *func, double *stackpos, uint8_t *color) {
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
        /*if (fabs(250*250*((x-xp)*(yp-ypp) - (xp-xpp)*(y-yp))) > 10) {
            printf("retrying\n");
            dt = dt/2;
            continue;
        }*/
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
}

void draw_function_slow(function *func, double *stackpos, uint8_t *color) {
    double x, y, xp, yp;
    eval_func(window_x0, &xp, &yp, func, stackpos);
    double minstep_x = (window_x1 - window_x0)/WIDTH;
    for (double t = window_x0+minstep_x; t < window_x1; t+=minstep_x) {
        eval_func(t, &x, &y, func, stackpos);
        LineAA((uint8_t*)data, WIDTH, HEIGHT, SCALE_XK(xp), SCALE_YK(yp), SCALE_XK(x), SCALE_YK(y), color);
        xp = x; yp = y;
    }
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

void redraw_all(gpointer data_pointer) {
    xscale = (1<<16)*1.0*WIDTH/(window_x1 - window_x0);
    yscale = (1<<16)*1.0*HEIGHT/(window_y1 - window_y0);
    clock_t t3, t4;
    t1 = clock();
    memset(data, 0, 3*WIDTH*HEIGHT);
    if ((window_x0 <= 0) && (0 < window_x1)) {
        int64_t xk = SCALE_XK(0) / (1<<16);
        for (uint32_t i=3*xk; i < 3*WIDTH*HEIGHT; i+=3*WIDTH) {
            data[i] = 128;
            data[i+1] = 128;
            data[i+2] = 128;
        }
    }
    if ((window_y0 < 0) && (0 <= window_y1)) {
        int64_t yk = SCALE_YK(0)/(1<<16);
        for (uint16_t i=0; i < 3*WIDTH; i++) {
            data[3*yk*WIDTH + i] = 128;
        }
    }
    
    t3 = clock();
    int64_t pt_x0, pt_x1, pt_y0, pt_y1;
    for (int i=0; i < n_expr; i++) {
        if (expression_list[i].flags & EXPRESSION_PLOTTABLE) {
            if ((expression_list[i].flags & EXPRESSION_FIXED) && (expression_list[i].value_type & TYPE_POINT)) {
                int len = (expression_list[i].value_type) >> 8;
                double *ptr = expression_list[i].value;
                for (int p=0; p < len; p+=2) {
                    pt_x0 = (SCALE_XK(ptr[p])>>16)-10;
                    pt_x1 = CLIP_WIDTH(pt_x0 + 20);
                    pt_x0 = CLIP_WIDTH(pt_x0);
                    pt_y0 = (SCALE_YK(ptr[p+1])>>16)-10;
                    pt_y1 = CLIP_HEIGHT(pt_y0 + 20);
                    pt_y0 = CLIP_HEIGHT(pt_y0);
                    g_print("point at (%f, %f), bounds %lld, %lld, %lld, %lld\n", ptr[p], ptr[p+1], pt_x0, pt_x1, pt_y0, pt_y1);
                    if (pt_x0 == pt_x1) continue;
                    for (int j=pt_y0; j < pt_y1; j++) {
                        for (int k=pt_x0; k < pt_x1; k++) memcpy(data + TOIDX(k, j), expression_list[i].color, 3);
                    }
                }
            } else {
                draw_function_slow(expression_list[i].func, stack+60, expression_list[i].color);
            }
        }
    }
    t4 = clock();
    if (data_pointer) {
        pixbuf = gdk_pixbuf_new_from_data(data, 0, 0, 8, WIDTH, HEIGHT, 3*WIDTH, NULL, NULL);
        gtk_image_set_from_pixbuf(data_pointer, pixbuf);
    }
    t2 = clock();
    g_print("Redraw took %luus, evaluation %luus, %d expressions\n", t2-t1, t4-t3, n_expr);
}

static gboolean scroll_callback (GtkWidget *event_box, GdkEventScroll *event, gpointer data_pointer) {
    //g_print ("Scrolled at %f, %f, direction %d, data_pointer: %p\n", event->x, event->y, event->direction, data_pointer);
    double center_x = ((int32_t)(event->x)<<16)/xscale + window_x0;
    double center_y = window_y1 - ((int32_t)(event->y)<<16)/yscale;
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
        redraw_all(data_pointer);
    }
    return TRUE;
}

static gboolean motion_callback(GtkWidget *event_box, GdkEventMotion *event, gpointer data_pointer) {
    if (click_state) {
        double dx = (event->x - click_x)/xscale*(1<<16);
        double dy = (click_y - event->y)/yscale*(1<<16);
        window_x0 -= dx; window_x1 -= dx; window_y0 -= dy; window_y1 -= dy;
        click_x = event->x;
        click_y = event->y;
        uint8_t color[4] = {255, 255, 255, 0};
        redraw_all(data_pointer);
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

void print_object(uint32_t type, double *pos) {
    if (type & TYPE_LIST) printf("[");
    uint32_t len = type >> 8;
    if ((type & TYPE_MASK) == TYPE_POINT) {
        for (int i=0; i < len; i += 2) {
            if (i+2 < len) printf("(%f, %f), ", pos[i], pos[i+1]);
            else printf("(%f, %f)", pos[i], pos[i+1]);
        }
    } else if ((type & TYPE_MASK) == TYPE_COLOR) {
        for (int i=0; i < len; i += 3) {
            if (i+2 < len) printf("rgb(%f, %f, %f), ", pos[i], pos[i+1], pos[i+2]);
            else printf("rgb(%f, %f, %f)", pos[i], pos[i+1], pos[i+2]);
        }
    } else {
        for (int i=0; i < len; i++) {
            if (i+1 < len) printf("%f, ", pos[i]);
            else printf("%f", pos[i]);
        }
    }
    if (type & TYPE_LIST) printf("]");
}

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

    pixbuf = gdk_pixbuf_new_from_data(data, 0, 0, 8, WIDTH, HEIGHT, 3*WIDTH, NULL, NULL);
    image = gtk_image_new_from_pixbuf(pixbuf);
    g_signal_connect(G_OBJECT(event_box), "button_press_event", G_CALLBACK(button_press_callback), image);
    g_signal_connect(G_OBJECT(event_box), "button_release_event", G_CALLBACK(button_release_callback), image);
    g_signal_connect(G_OBJECT(event_box), "scroll_event", G_CALLBACK(scroll_callback), image);
    g_signal_connect(G_OBJECT(event_box), "motion-notify-event", G_CALLBACK(motion_callback), image);
    gtk_container_add(GTK_CONTAINER(event_box), image);

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
    parse_file(function_list, stack, variable_list, stringbuf, expression_list, &n_func, &n_var, &n_expr);
    printf("Parsing completed. %d function blocks, %d variables, %d expressions\n", n_func, n_var, n_expr);

    for (int i=0; i < n_expr; i++) {
        if ((i%7+1)&0x01) expression_list[i].color[0] = 255;
        else expression_list[i].color[0] = 0;
        if ((i%7+1)&0x02) expression_list[i].color[1] = 255;
        else expression_list[i].color[1] = 0;
        if ((i%7+1)&0x04) expression_list[i].color[2] = 255;
        else expression_list[i].color[2] = 0;
    }

    uint32_t type;
    t1 = clock();
    for (int i=0; i < 1; i++) 
        type = expression_list[5].func->oper(expression_list[5].func, stack+60);
    t2 = clock();
    printf("result is ");
    print_object(type, stack+60);
    printf(", took %luus\n", t2 - t1);
    print_object(expression_list[5].value_type, expression_list[5].value); printf("\n");

    type = expression_list[6].func->oper(expression_list[6].func, stack+60);
    printf("result is ");
    print_object(type, stack+60);
    printf("\n");

    stack[60] = 37;
    variable_list[0].pointer = stack+60;
    variable_list[0].type = 1<<8;
    printf("result is ");
    print_object(expression_list[7].func->oper(expression_list[7].func, stack+61), stack+61);
    printf(", expression flags %02x, value type %08x\n", expression_list[7].flags, expression_list[7].value_type);

    //return 0;

    /*t1 = clock();

    for (int i=0; i < 500; i++) {
        array[i] = i;
    }
    variable_list[0].pointer = array;
    variable_list[0].type = 1000<<8;
    printf("variable_list %p, %p, array %p, %p\n", &(variable_list[0].type), &(variable_list[0].pointer), array, &(array[1]));
    function_list[0].oper(function_list, stack+10);
    t2 = clock();
    for (int i=0; i < 10; i++)
        printf("Result is %f, %luus\n", stack[10+i], t2-t1);
    return 0;*/

    GtkApplication *app;
    int status;
    redraw_all(NULL);

    /*t = t + dt;
    draw_line(SCALE_X(xpp), SCALE_Y(ypp), SCALE_X(xp), SCALE_Y(yp));
    while (t < 1) {
        newdt = 0;
        eval_func(t+dt, &x, &y);
        est_radius(x, y, xp, yp, xpp, ypp, &radius, &speed);
        newdt = fmin(RADIUS_SCALE * radius, MAX_ARC_LENGTH * dt / speed);
        printf("t: %f -> (%f, %f), radius: %f, arc: %f, newdt: %f, dt: %f, det: %f\n", t, x, y, radius, speed, newdt, dt, (x-xp)*(yp-ypp) - (xp-xpp)*(y-yp));
        if (newdt*1.5 < dt) {
            printf("  retrying\n");
            // Retry with a smaller step if the curvature would be too high with the current step.
            dt = newdt;
            continue;
        }
        //if (newdt > dt) dt = newdt;
        draw_line(SCALE_X(xp), SCALE_Y(yp), SCALE_X(x), SCALE_Y(y));
        t = t + dt;
        dt = newdt;
        xpp = xp; ypp = yp; xp = x; yp = y;
    }*/


    //draw_line(200, 200, 200, 250);
    //draw_line(200, 250, 200, 300);
    /*for (int i=0; i < 12; i++) {
        double angle = i * 3.1415926 / 6;
        draw_line(200, 200, 200 + (int)(50*cos(angle)), 200 + (int)(50*sin(angle)));
        draw_line(200 + (int)(50*cos(angle)), 200 + (int)(50*sin(angle)), 200 + (int)(100*cos(angle)), 200 + (int)(100*sin(angle)));
    }*/

    app = gtk_application_new("org.gtk.example", G_APPLICATION_FLAGS_NONE);
    g_signal_connect(app, "activate", G_CALLBACK (activate), NULL);
    status = g_application_run(G_APPLICATION (app), argc, argv);
    g_object_unref (app);

    return status;
}
