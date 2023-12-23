#include <stdint.h>
#include <gtk/gtk.h>

#ifndef PARSE_H
#define PARSE_H

#define VARIABLE_IN_SCOPE 0x01
#define VARIABLE_ARGUMENT 0x02
#define VARIABLE_FUNCTION 0x04
#define VARIABLE_SPECIAL 0x08
#define VARIABLE_XLIKE 0x10
#define VARIABLE_YLIKE 0x20
#define VARIABLE_XYLIKE 0x30
#define VARIABLE_ACTION 0x40
#define VARIABLE_INTERVAL 0x80

#define EXPRESSION_PLOTTABLE 0x01
#define EXPRESSION_FIXED 0x02
#define EXPRESSION_ACTION 0x04

#define PARSE_COMMA 0x01
#define PARSE_NEWVAR 0x02
#define PARSE_ACTION 0x04

#define FIND_NEAREST_UNSUPPORTED 0x01
#define FIND_NEAREST_NOCONV 0x02

/*
 * Data structures used for this program
 *  * Function list that stores function blocks
 *  * Stack that stores constants, variables, and temporary results
 *  * Variable list that stores variable blocks
 *  * Expression list that stores expression blocks
 *  * List of groups of expression pointers that depend on a particular expression
 */

typedef struct function_s {
    uint32_t (*oper)(void*, double*);       // Pointer to the function this function block implements. May not be null
    uint32_t (*inter)(void*, double*, double*);
    void *value;                            // Points to the location of a variable. This is used when oper is func_value. May be null otherwise
                                            // Also used to point to a function definition when the operation is another user-defined function
    uint32_t value_type;
    struct function_s *next_arg;            // Points to the parent function's next argument, null otherwise
    struct function_s *first_arg;           // Points to the first child
} function;

typedef struct variable_s {
    char *name;
    uint8_t flags;
    double *pointer;
    uint32_t type;
    double *new_pointer;
    uint32_t new_type;
} variable;

typedef struct expression_s {
    struct variable_s *var;                 // Variable to which the expression is assigned, may be NULL
    struct function_s *func;                // Root operation of the expression
    struct expression_s **dependencies;     // List of variables that this expression depends upon
    uint8_t num_dependencies;               // Number of variables this expression depends upon
    uint8_t num_nonfixed_dependencies;      // Number of non-fixed variables this expression depends upon
    uint8_t num_dependents;                 // Number of expressions that depend on this expression.
    int expr_begin;                         // Position on the line where the definition of the expression begins
    struct expression_s *next_expr;         // Next expression to be evaluated. Assigned during sort.
    uint8_t flags;                          // Flags for the expression
    uint8_t color[4];                       // Color in RGBA, 0-255
    struct function_s *color_pointer;       // Pointer to a function block that produces the color
    double *value;                          // Pointer to the value of the expression. Used when plotting
    uint32_t value_type;                    // Type of the last value of the expression. Used when plotting
    char *def;                              // Pointer to char array containing definition. Freed after parsing
} expression;

typedef struct file_data_s {
    struct expression_s *expression_list;
    struct variable_s *variable_list;
    struct function_s *function_list;
    double *stack;
    double *lstack;
    uint32_t n_expr;
    uint32_t n_var;
    uint32_t n_func;
    uint32_t n_stack;
    GtkWidget *drawing_area;
    struct expression_s **deptable;
    cairo_surface_t *surface;
    cairo_surface_t *overlay;
    uint8_t use_overlay;
    struct expression_s *click_expr;
} file_data;

// Struct that contains information about a particular operator
typedef struct oper_data_s {
    uint32_t (*oper)(void*, double*);           // Pointer to the operator function
    char *name;                                 // Name of the operator
    uint32_t (*inter)(void*, double*, double*); // Pointer to the operator's interval function
} oper_data;

typedef struct integration_result_s {
    struct function_s *func;                // Function being integrated
    double lbv;                             // Lower integration bound
    double ubv;                             // Upper integration bound
    double result;                          // Result of the integration
    uint8_t negate;                         // 1 if the result of the integration should be negated
} integration_result;

double parse_double(char *string);

function new_function(uint32_t (*oper)(void*, double*), function *next_arg, function *first_arg);
function new_value(void *value, uint32_t value_type, function *next_arg);
variable new_variable(char *name, int type, uint8_t flags, double *pointer);

void print_object(uint32_t type, double *pos);

void load_file(char *fname, file_data *fd);
void evaluate_from(file_data *fd, expression *top_expr);
int parse_latex_rec(char *latex, int end, function *function_list, double *stack, variable *variable_list, char *stringbuf, int *stack_size, int *var_size, int *string_size, uint8_t *flags);
void parse_latex(char *latex, function *function_list, double *stack, variable *variable_list, char *stringbuf, int *stack_size, int *var_size, int *string_size);
expression* parse_file(file_data *fd, char *stringbuf);

#endif
