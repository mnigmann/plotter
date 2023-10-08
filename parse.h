#include <stdint.h>

#ifndef PARSE_H
#define PARSE_H

#define VARIABLE_IN_SCOPE 0x01
#define VARIABLE_ARGUMENT 0x02
#define VARIABLE_FUNCTION 0x04
#define VARIABLE_SPECIAL 0x08
#define VARIABLE_XLIKE 0x10
#define VARIABLE_YLIKE 0x20
#define VARIABLE_XYLIKE 0x30

#define EXPRESSION_PLOTTABLE 0x01
#define EXPRESSION_FIXED 0x02

#define PARSE_COMMA 0x01
#define PARSE_NEWVAR 0x02

/*
 * Data structures used for this program
 *  * Function list that stores function blocks
 *  * Stack that stores constants, variables, and temporary results
 *  * Variable list that stores variable blocks
 *  * Expression list that stores expression blocks
 *  * List of groups of expression pointers that depend on a particular expression
 */

typedef struct function_s {
    uint32_t (*oper)(void*, double*);   // Pointer to the function this function block implements
    void *value;                        // Points to the location of a variable. This is used when oper is null. May be null otherwise
                                        // Also used to point to a function definition when the operation is another user-defined function
    uint32_t value_type;
    struct function_s *next_arg;  // Points to the parent function's next argument, or the parent function if this is the last argument. Null in the root function
    struct function_s *first_arg; // Points to the first child
} function;

typedef struct variable_s {
    char *name;
    uint32_t type;
    uint8_t flags;
    double *pointer;
} variable;

typedef struct expression_s {
    struct variable_s *var;                 // Variable to which the expression is assigned, may be NULL
    struct function_s *func;                // Root operation of the expression
    struct expression_s **dependencies;     // List of variables that this expression depends upon
    uint8_t num_dependencies;               // Number of variables this expression depends upon
    uint8_t num_dependents;
    int expr_begin;                         // Position on the line where the definition of the expression begins
    struct expression_s *next_expr;         // Next expression to be evaluated. Assigned during sort.
    uint8_t flags;                          // Flags for the expression
    uint8_t color[4];
    double *value;
    uint32_t value_type;
} expression;

function new_function(uint32_t (*oper)(void*, double*), function *next_arg, function *first_arg);
function new_value(void *value, uint32_t value_type, function *next_arg);
variable new_variable(char *name, int type, uint8_t flags, double *pointer);

void print_object(uint32_t type, double *pos);

void evaluate_from(expression *expression_list, int n_expr, expression *top_expr, double *stack, int *stack_size_ptr);
int parse_latex_rec(char *latex, int end, function *function_list, double *stack, variable *variable_list, char *stringbuf, int *stack_size, int *var_size, int *string_size, uint8_t *flags);
void parse_latex(char *latex, function *function_list, double *stack, variable *variable_list, char *stringbuf, int *stack_size, int *var_size, int *string_size);
int parse_file(char *fname, function *function_list, double *stack, variable *variable_list, char *stringbuf, expression *expression_list, uint32_t *n_func, uint32_t *n_var, uint32_t *n_expr);

#endif
