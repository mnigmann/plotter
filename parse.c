#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "parse.h"
#include "functions.h"
#include "linalg_functions.h"
#include <math.h>

#define N_FUNCTIONS 16
const char *function_names[N_FUNCTIONS] = {
    "arctan",       "cos",          "sin",          "mod",          "floor",        "polygon",      "total",        "distance",
    "rgb",          "max",          "abs",          "sort",         "join",         "length",       "solve",        "eigvals"
};
uint32_t (*function_pointers[N_FUNCTIONS])(void*, double*) = {
    func_arctan,    func_cosine,    func_sine,      func_mod,       func_floor,     func_polygon,   func_total,     func_distance,
    func_rgb,       func_max,       func_abs,       func_sort,      func_join,      func_length,    func_solve,     func_eigvals
};

int extract_braces(char* latex, int start) {
    int end = strlen(latex);
    int braces = 0;
    for (int i=start; i < end; i++) {
        if (latex[i] == '{') braces++;
        if (latex[i] == '}') braces--;
        if (braces == 0) return i;
    }
    return -1;
}

int extract_brackets(char *latex, int start) {
    int end = strlen(latex);
    int parentheses = 0;
    for (int i=start; i+6 < end; i++) {
        if (strncmp(latex+i, "\\left[", 6) == 0) {
            parentheses++;
            i += 5;
        } else if (strncmp(latex+i, "\\right]", 7) == 0) {
            parentheses--;
            i += 6;
        }
        if (parentheses == 0) return i;
    }
    return -1;
}

int extract_parenthetical(char *latex, int start) {
    int end = strlen(latex);
    int parentheses = 0;
    for (int i=start; i+6 < end; i++) {
        if (strncmp(latex+i, "\\left(", 6) == 0) {
            parentheses++;
            i += 5;
        } else if (strncmp(latex+i, "\\right)", 7) == 0) {
            parentheses--;
            i += 6;
        }
        if (parentheses == 0) return i;
    }
    return -1;
}

int extract_double(char *latex, int start, double *result) {
    int end = strlen(latex);
    uint8_t is_constant = 1;
    double pow = 0.1;
    double value = 0;
    for (int i=start; i < end; i++) {
        if ((i == start) && (latex[i] == '-')) is_constant |= 0x02;
        else if (latex[i] == '.') is_constant |= 0x04;
        else if (('0' <= latex[i]) && (latex[i] <= '9')) {
            if (is_constant & 0x04) {
                value += pow*(latex[i] - '0');
                pow = pow/10;
            } else value = 10*value + (latex[i] - '0');
        } else {
            *result = (is_constant & 0x02 ? -value : value);
            return i-1;
        }
    }
    *result = (is_constant & 0x02 ? -value : value);
    return end-1;
}

double parse_double(char *string) {
    double res = NAN;
    extract_double(string, 0, &res);
    return res;
}

function new_function(uint32_t (*oper)(void*, double*), function *next_arg, function *first_arg) {
    function block;
    block.oper = oper;
    block.value = NULL;
    block.value_type = 0;
    block.next_arg = next_arg;
    block.first_arg = first_arg;
    return block;
}

function new_value(void *value, uint32_t value_type, function *next_arg) {
    function block;
    block.oper = func_value;
    block.value = value;
    block.value_type = value_type;
    block.next_arg = next_arg;
    block.first_arg = NULL;
    return block;
}

variable new_variable(char *name, int type, uint8_t flags, double *pointer) {
    variable var;
    var.name = name;
    var.type = type;
    var.flags = flags;
    var.pointer = pointer;
    var.new_pointer = NULL;
    var.new_type = 0;
    return var;
}

expression new_expression(variable *var, function *func, expression **dependencies, uint8_t num_dependencies, int expr_begin) {
    expression expr;
    expr.var = var;
    expr.func = func;
    expr.dependencies = dependencies;
    expr.num_dependencies = num_dependencies;
    expr.num_dependents = 0;
    expr.flags = 0;
    expr.expr_begin = expr_begin;
    expr.next_expr = NULL;
    expr.value = NULL;
    expr.value_type = 0;
    return expr;
}

void shift_blocks(function *function_list, uint32_t start, uint32_t num) {
    function *start_addr = function_list + start;
    function *end_addr = function_list + (start+num);
    for (int j=start+num-1; j >= (int)start; j--) {
        function_list[j+1] = function_list[j];
        if ((start_addr <= function_list[j+1].next_arg) && (function_list[j+1].next_arg < end_addr)) function_list[j+1].next_arg++;
        if ((start_addr <= function_list[j+1].first_arg) && (function_list[j+1].first_arg < end_addr)) function_list[j+1].first_arg++;
    }
}

int insert_product_term(function *function_list, uint32_t last_pos, uint32_t func_pos) {
    if (func_pos) {
        if (last_pos == 0) {
            // If last_pos is zero but func_pos is not, then one product term has already been
            // inserted. Thus, we must insert a multiplication block
            shift_blocks(function_list, last_pos, func_pos-last_pos);
            function_list[last_pos] = new_function(func_multiply, function_list[last_pos+1].next_arg, function_list+last_pos+1);
            last_pos++;
            func_pos++;
        }
        function_list[last_pos].next_arg = function_list + func_pos;
    }
    return func_pos;
}

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

void print_object_short(uint32_t type, double *pos) {
    if ((type & TYPE_LIST) && ((type>>8) > 100)) {
        uint32_t len = (type>>8);
        if ((type & TYPE_MASK) == TYPE_POINT) len = len / 2;
        if ((type & TYPE_MASK) == TYPE_COLOR) len = len / 3;
        printf("[ %d total elements, type %08x ]\n", len, type);
        return;
    }
    if (type & TYPE_LIST) printf("[");
    uint32_t len = type >> 8;
    if ((type & TYPE_MASK) == TYPE_POINT) {
        for (int i=0; i < len; i += 2) {
            if (i+2 < len) printf("(%0.12f, %0.12f), ", pos[i], pos[i+1]);
            else printf("(%0.12f, %0.12f)", pos[i], pos[i+1]);
        }
    } else if ((type & TYPE_MASK) == TYPE_COLOR) {
        for (int i=0; i < len; i += 3) {
            if (i+2 < len) printf("rgb(%f, %f, %f), ", pos[i], pos[i+1], pos[i+2]);
            else printf("rgb(%f, %f, %f)", pos[i], pos[i+1], pos[i+2]);
        }
    } else {
        for (int i=0; i < len; i++) {
            if (i+1 < len) printf("%0.12f, ", pos[i]);
            else printf("%0.12f", pos[i]);
        }
    }
    if (type & TYPE_LIST) printf("]");
}

/*
int extract_variable(char *start, int length, function *function_list, variable *variable_list) {
    int varindex = -1;
    uint8_t flags;
    for (int j=0; variable_list[j].name; j++) {
        uint8_t flags = variable_list[j].flags;
        if ((flags & VARIABLE_IN_SCOPE) && (strncmp(variable_list[j].name, start, length) == 0)) {
            varindex = j;
            if (flags & VARIABLE_ARGUMENT) break; // If an argument is in scope, prefer it
        }
    }
    if (varindex == -1) {
        printf("ERROR: variable %.*s not found!\n", length, start);
        exit(EXIT_FAILURE);
    }
    if (variable_list[varindex].flags & VARIABLE_FUNCTION) {
        // Function found
        // We cannot assume that the function has been defined, only declared. The function
        // definitions must be correctly connected later. For now, we connect to the variable
        // declaration
        int arg1_end = extract_parenthetical(start, length);
        printf("multiply %.*s at %d with argument %.*s\n", cmd_len+1, latex+cmd_start-1, varindex, arg1_end-cmd_end+1, latex+cmd_end);
        func_pos = insert_product_term(function_list, last_pos, func_pos);
        function_list[func_pos] = new_function(func_user_defined, NULL, function_list+func_pos+1);
        function_list[func_pos].value = variable_list+varindex;
        last_pos = func_pos;
        func_pos++;
        func_pos += PARSE_LATEX_REC(latex+cmd_end+6, arg1_end-cmd_end-12, function_list+func_pos);
        i = arg1_end;
    } else {
        // Variable found
        printf("multiply %.*s at %d, func_pos %d\n", cmd_len+1, latex+cmd_start-1, varindex, func_pos);
        func_pos = insert_product_term(function_list, last_pos, func_pos);
        printf("func_pos %d\n", func_pos);
        function_list[func_pos] = new_value(variable_list+varindex, 0x40, NULL);
        last_pos = func_pos;
        func_pos++;
        i = cmd_end-1;
    }
}*/

int get_next_match(char *latex, int start, int end, char target) {
    int n_braces = 0;
    int n_leftright = 0;
    for (int i=start; i < end; i++) {
        if (latex[i] == '{') n_braces++;
        else if (latex[i] == '}') n_braces--;
        else if (strncmp(latex+i, "\\left", 5) == 0) {
            n_leftright++;
            i += 4;
        } else if (strncmp(latex+i, "\\right", 6) == 0) {
            n_leftright--;
            i += 5;
        } else if ((latex[i] == target) && (n_leftright == 0) && (n_braces == 0)) {
            return i;
        }
    }
    return -1;
}

void evaluate_from(file_data *fd, expression *top_expr) {
    // Print the expression table
    //for (int i=0; i < n_expr; i++)
    //    printf("expression pointer at %p, function %p, variable %p, flags %02x, begin %d, variable %p, type %08x\n", expression_list+i, expression_list[i].func, expression_list[i].var, expression_list[i].flags, expression_list[i].expr_begin, expression_list[i].value, expression_list[i].value_type);


    uint32_t type;
    expression *expr = top_expr;
    expression *expression_list = fd->expression_list;
    double *stack = (fd->stack) + (fd->n_stack);
    uint32_t n_expr = fd->n_expr;
    int stack_size = -1;
    //printf("stack_size is %d with value %p, stack %p\n", stack_size, expr->value, stack);
    double *ptr;
    while (expr) {
        if ((expr->func) && (expr->var) && !(expr->var->flags & (VARIABLE_FUNCTION | VARIABLE_ACTION)) && ((expr->flags & EXPRESSION_FIXED) || !(expr->flags & EXPRESSION_PLOTTABLE))) {
#ifdef DEBUG_EVAL
            printf("evaluating variable %s, expression %p (%03ld), variable block %p, function block %p, old pointer %p, stack %p\n", expr->var->name, expr, expr-expression_list+1, expr->var, expr->func, expr->value, stack);
#endif
            type = (expr->func->oper(expr->func, stack));
            ptr = NULL;
            // If the size changes, reallocate. Otherwise, use the same memory
            if (((expr->value_type) >> 8) != (type >> 8)) {
                if ((expr->value_type) >> 8) ptr = realloc(expr->value, (type>>8)*sizeof(double));
                else ptr = malloc((type>>8)*sizeof(double));
            } else ptr = expr->value;
            // Special case for zero-length lists
            // TODO: actions
            if (type>>8) memcpy(ptr, stack, (type>>8)*sizeof(double));
            else if (!ptr) ptr = stack;
#ifdef DEBUG_EVAL
            printf("    result "); print_object_short(type, ptr);
            printf(" stored to %p, has type %08x\n", ptr, type);
#endif
            if (type & TYPE_POINT) expr->flags |= EXPRESSION_PLOTTABLE | EXPRESSION_FIXED;
            expr->var->pointer = ptr;
            expr->var->type = type;
            expr->value = ptr;
            expr->value_type = type;
        } //else if (expr->value) stack_size += ((expr->value_type) >> 8);
        expr = expr->next_expr;
    }

    // Evaluate all expressions that are not definitions or actions and are not dependent on x or y
    for (expr=expression_list; expr < expression_list+n_expr; expr++) {
        if ((expr->flags & EXPRESSION_FIXED) && !(expr->var) && !(expr->flags & EXPRESSION_ACTION)) {
#ifdef DEBUG_EVAL
            printf("evaluating expression %p (%03ld)\n", expr, expr-expression_list+1);
#endif
            type = (expr->func->oper(expr->func, stack));
            // If the size changes, reallocate. Otherwise, use the same memory
            if (((expr->value_type) >> 8) != (type >> 8))
                ptr = realloc(expr->value, (type>>8)*sizeof(double));
            else ptr = expr->value;
            memcpy(ptr, stack, (type>>8)*sizeof(double));
            if (((type & TYPE_MASK) == TYPE_POINT) || ((type & TYPE_MASK) == TYPE_POLYGON)) expr->flags |= EXPRESSION_PLOTTABLE;
#ifdef DEBUG_EVAL
            printf("    result "); print_object_short(type, ptr);
            printf(" stored to %p, has type %08x\n", ptr, type);
#endif
            expr->value = ptr;
            expr->value_type = type;
        }
        if (!(expr->func) && (expr->var)) {
#ifdef DEBUG_EVAL
            printf("expression %p is variable %s with value ", expr, expr->var->name); print_object_short(expr->var->type, expr->var->pointer); printf("\n");
#endif
            expr->value = expr->var->pointer;
            expr->value_type = expr->var->type;
        }
    }
}

/*
L_{endpoints}=N_{c}\left[A_{pindex}\left(
    N_{f},
    G\left(N_{s0}\left[A_{csnnn}\right],
        N_{s1}\left[A_{csnnn}\right],
        \left[1...\operatorname{length}\left(A_{csnnn}\right)\right]-A_{csnn}\left[A_{csnnn}\right]-1,
        N_{st}\left[A_{csnnn}\right]
    \right)
\right)\right]
 */

/* LaTeX operators:
 *  addition (+), subtraction (-), multiplication (\cdot), exponentiation (^), subscript(_)
 *  Exponentiation always uses braces on the exponent, maybe parentheses on the base
 
 *  A\left[B\right] --> index(A, B) --> double or list
 *  a...b           --> range(a, b) --> list
 
 * Things that can be mutliplied together:
 *  Exponential expressions (a^{b})
 *  Parenthetical expressions
 *  Variables
 *  Indexed variables (A\left[B\right])
 *  Lists (\cdot\left[A\right])
 *  
 *
 *
 * Creating callable functions
 *  A callable function must have a specific format. Consider a function with n arguments.
 *  The outermost block of the function is stored in the lowest position. The next n blocks are value blocks
 *  into which the argument types, values, and sizes are stored. The remaining blocks are function code.
 * 
 * Global variables
 *  A special table is needed to store the string names of global variables and their addresses in the memory.
 *  Each value block that references a variable instead references an index in this table, allowing for
 *  variables that change size.
 *
 * Function table
 *  Another table stores the string names of functions and their starting function block.
 *
 *
 * Process
 *  Extract variables
 *  Determine independent and dependent variables
 *  Check if the number of free variables is less than or equal to one.
 *  Replace that variable with x
 *
 *  In the case of a regression (~ operator), we need to define additional free variables that are bound to the
 *  regression. 
 *
 *
 *
 * We will need a graph storing the relationships between variables. If there is a circular definition, we
 * need to raise an error. This goes for functions as well.
 * 
 * 
 * FEATURES WE NEED TO IMPLEMENT:
 *  * Proper support for \pi, \tau, \alpha, \beta, \theta
 *  * Constants \pi, \tau, e
 *  * Support for points and point operations
 *  * Store the value of expressions independent of x instead of re-evaluating them when we re-draw
 *  * Actions
 *  * Conditionals
 *  * Correctly identify what must be updated when something changes.
 *  * Evaluate expressions with 0 dependencies ahead of time and remove them from the update list.
 *  * for statements of the form \left[expr1\operatorname{for}var=expr2\right]
 *      * The variable(s) must be temporarily brought into scope when parsing expr1 and removed immediately after
 *
 * TYPES
 *  * double (1)
 *      * Any function
 *  * point (2)
 *      * Multiply by scalar
 *      * Add or subtract another point
 *      * Extract one of the components
 *      * Compute distance
 *  * color (3)
 *      * No operations are supported on colors
 *  * boolean (1)
 *      * May be used as an index or as the input to a conditional
 *  * list of any of the above
 *
 * PLOTTABILITY
 * An expression can be plotted if:
 *  * the expression is an expression and ultimately only depends on x or y
 *  * the expression is an equation and ultimately depends on x and y
 *  * the expression is an equation defining r in terms of \theta
 *  * the expression is a variable definition and has exactly one free variable
 *  * the expression is a function definition and the function has exactly one argument
 *  * the expression is a parametric expression with exactly one free variable
 *  * the expression is a point or list of points
 */

#define PARSE_LATEX_REC(_latex, _end, _function_list) parse_latex_rec(_latex, _end, _function_list, stack, variable_list, stringbuf, stack_size, var_size, string_size, &flags)

int parse_latex_rec(char *latex, int end, function *function_list, double *stack, variable *variable_list, char *stringbuf, int *stack_size, int *var_size, int *string_size, uint8_t *result_flags) {
    printf("Parsing %.*s\n", end, latex);
    
    if (end == 0) return 0;

    // Check if the expression is an ellipsis
    if (end == 3 && strncmp(latex, "...", 3) == 0) {
        function_list[0] = new_function(func_ellipsis, NULL, NULL);
        return 1;
    }
    // Check if the expression is a constant
    double value = 0;
    int double_len = extract_double(latex, 0, &value);
    if (double_len+1 == end) {
        stack[*stack_size] = value;
        function_list[0] = new_value(stack+(*stack_size), (1<<8), NULL);
        *stack_size += 1;
        printf("Found constant %f\n", value);
        return 1;
    }

    int last_term = 0;
    uint8_t subtract = 0;
    uint8_t n_leftright = 0;
    uint8_t n_braces = 0;
    uint8_t n_terms = 0;
    int func_pos = 0;
    int last_pos = 0;
    uint8_t flags;
    int subscript, superscript, start;

    // Decompose commas
    last_term = 0;
    for (int i=0; i < end; i++) {
        i = get_next_match(latex, i, end, ',');
        if (n_terms && (i == -1)) i = end;
        if (i >= 0) {
            printf("comma term %.*s\n", i - last_term, latex+last_term);
            
            if (n_terms > 0) {
                // For assignments, we chain by value and not by next_arg
                printf("Chaining comma term with flags %02x\n", flags);
                function_list[last_pos].next_arg = function_list + func_pos;
            }
            last_pos = func_pos;
            flags = 0;
            func_pos += PARSE_LATEX_REC(latex+last_term, i - last_term, function_list+func_pos);
            last_term = i+1;
            n_terms++;
        } else break;
    }
    if (n_terms) {
        *result_flags |= PARSE_COMMA;
        printf("parsing %.*s done, func_pos %d\n", end, latex, func_pos);
        return func_pos;
    }

    // Decompose equality (=)
    n_braces = 0;
    n_leftright = 0;
    int index_equals = get_next_match(latex, 0, end, '=');
    if (index_equals >= 0) {
        function_list[0] = new_function(func_equals, NULL, function_list+1);
        func_pos++;
        func_pos += PARSE_LATEX_REC(latex, index_equals, function_list+func_pos);
        function_list[1].next_arg = function_list + func_pos;
        index_equals++;
        if (latex[index_equals] == ' ') index_equals++;
        func_pos += PARSE_LATEX_REC(latex+index_equals, end-index_equals, function_list+func_pos);
        return func_pos;
    }
    n_braces = 0;
    n_leftright = 0;
    index_equals = get_next_match(latex, 0, end, '>');
    if (index_equals >= 0) {
        function_list[0] = new_function(func_greater, NULL, function_list+1);
        func_pos++;
        func_pos += PARSE_LATEX_REC(latex, index_equals, function_list+func_pos);
        function_list[1].next_arg = function_list + func_pos;
        index_equals++;
        if (latex[index_equals] == ' ') index_equals++;
        func_pos += PARSE_LATEX_REC(latex+index_equals, end-index_equals, function_list+func_pos);
        return func_pos;
    }

    // Decompose actions
    for (int i=0; i < end; i++) {
        if (latex[i] == '{') n_braces++;
        else if (latex[i] == '}') n_braces--;
        else if (strncmp(latex+i, "\\left", 5) == 0) {
            n_leftright++;
            i += 4;
        } else if (strncmp(latex+i, "\\right", 6) == 0) {
            n_leftright--;
            i += 5;
        } else if ((strncmp(latex+i, "\\to", 3) == 0) && (n_leftright == 0) && (n_braces == 0)) {
            function_list[0] = new_function(func_assign, NULL, function_list+1);
            func_pos++;
            func_pos += PARSE_LATEX_REC(latex, i, function_list+func_pos);
            function_list[1].next_arg = function_list + func_pos;
            i += 3;
            if (latex[i] == ' ') i++;
            func_pos += PARSE_LATEX_REC(latex+i, end-i, function_list+func_pos);
            *result_flags |= PARSE_ACTION;
            return func_pos;
        }
    }
        

    // Decompose addition
    n_terms = 0;
    for (int i=0; i < end; i++) {
        if (latex[i] == '{') n_braces++;
        else if (latex[i] == '}') n_braces--;
        else if (strncmp(latex+i, "\\left", 5) == 0) {
            n_leftright++;
            i += 4;
        } else if (strncmp(latex+i, "\\right", 6) == 0) {
            n_leftright--;
            i += 5;
        } else if (((latex[i] == '+') || (latex[i] == '-')) && (n_leftright == 0) && (n_braces == 0)) {
            if (n_terms == 0) {
                function_list[0] = new_function(func_add, NULL, NULL);
                last_pos = 0;
                func_pos++;
            }
            if (subtract) printf("subtract %.*s\n", i - last_term, latex+last_term);
            else printf("add %.*s\n", i - last_term, latex+last_term);
            
            if (last_pos == 0) function_list[last_pos].first_arg = function_list + func_pos;
            else function_list[last_pos].next_arg = function_list + func_pos;
            last_pos = func_pos;
            
            func_pos += PARSE_LATEX_REC(latex+last_term, i - last_term, function_list+func_pos);
            if (subtract) function_list[last_pos].value_type ^= 0x80;
            subtract = (latex[i] == '-' ? 1 : 0);
            last_term = i+1;
            n_terms++;
        }
    }
    if (n_terms) {
        // If the expression consists of several things added together, we recurse on all of them. Otherwise, we continue parsing.
        if (subtract) printf("subtract (final) %.*s\n", end - last_term, latex+last_term);
        else printf("add (final) %.*s\n", end - last_term, latex+last_term);
        
        function_list[last_pos].next_arg = function_list + func_pos;
        last_pos = func_pos;
            
        func_pos += PARSE_LATEX_REC(latex+last_term, end - last_term, function_list+func_pos);
        if (subtract) function_list[last_pos].value_type ^= 0x80;
        printf("parsing %.*s done, func_pos %d\n", end, latex, func_pos);
        return func_pos;
    }
    
    // Decompose multiplication
    // At this point, we can assume that the expression does not consist of several things
    // added together. We can also assume that it does not begin with a negative sign
    // 
    // A product term may begin with either:
    //  * Numbers
    //  * \cdot
    //  * a variable name (including \pi, \tau, \alpha, \beta, and \theta), possibly with a subscript
    //  * a parenthetical expression of the form \left(:::\right)
    //  * a conditional expression of the form \left\{:::\right\}
    //  * a list (only for the first product term)
    //  * a point
    //  * an expression of the form \oper{x} or \frac{x}{y} or \operatorname{oper}\left(:::\right) or F\left(:::\right)
    // A product term may be followed by any number of exponents and indices in any order
    // If a second term is found, then the function list must be shifted. This eliminates possible extraneous multiplications

    // A_{x}^{y} --> index=0, subscript = 4, superscript = 8
    printf("Decomposing multiplication\n");
    uint8_t in_command = 0;
    int idx_end = 0;
    int cmd_end, cmd_len, cmd_start, arg1_end, arg1_start, arg2_end, arg2_start;
    func_pos = 0;
    uint32_t (*oper)(void*, double*);
    for (int i=0; i < end; i++) {
        // Check if the previous term has an index. If it does, we must apply the index by moving
        // the term code down one function block and inserting the index block above it.
        if ((func_pos > 0) && (i+5 < end) && (strncmp(latex+i, "\\left[", 6) == 0)) {
            // Variable has an index
            idx_end = extract_brackets(latex, i);
            printf("Index found at position %d - %d, func_pos %d, last_pos %d\n", i, idx_end, func_pos, last_pos);
            flags = 0;
            // Move the previous term down one block
            //for (int j=func_pos-1; j >= last_pos; j--) function_list[j+1] = function_list[j];
            shift_blocks(function_list, last_pos, func_pos-last_pos);
            func_pos++;
            // Parse the index
            uint32_t parsed_length = PARSE_LATEX_REC(latex+i+6, idx_end - i - 12, function_list+func_pos);
            // Insert index block and chain the term to the index
            function_list[last_pos] = new_function(func_index, function_list[last_pos+1].next_arg, function_list+last_pos+1);
            function_list[last_pos+1].next_arg = function_list+func_pos;
            if (flags & PARSE_COMMA) {
                // The index is a comma-separated list of terms. We need to add in a list block
                // This is done by shifting the index chain down
                //for (int j=func_pos+parsed_length-1; j >= func_pos; j--) function_list[j+1] = function_list[j];
                shift_blocks(function_list, func_pos, parsed_length);
                function_list[func_pos] = new_function(func_list, NULL, function_list+func_pos+1);
                func_pos += parsed_length+1;
            } else {
                // The index is a single expression and is connected directly to the index block
                func_pos += parsed_length;
            }
            i = idx_end;
        }
        // Check if the previous term has a superscript
        else if ((func_pos > 0) && (latex[i] == '^')) {
            superscript = extract_braces(latex, i+1);
            printf("superscript found %.*s\n", superscript-i-2, latex+i+2);
            shift_blocks(function_list, last_pos, func_pos-last_pos);
            func_pos++;
            function_list[last_pos] = new_function(func_exponentiate, function_list[last_pos+1].next_arg, function_list+last_pos+1);
            function_list[last_pos+1].next_arg = function_list+func_pos;
            func_pos += PARSE_LATEX_REC(latex+2+i, superscript - i - 2, function_list+func_pos);
            i = superscript;
        }
        // Check if the previous term has a coordinate (.x or .y)
        else if ((func_pos > 0) && (latex[i] == '.')) {
            shift_blocks(function_list, last_pos, func_pos-last_pos);
            func_pos++;
            if (latex[i+1] == 'x') {
                function_list[last_pos] = new_function(func_extract_x, function_list[last_pos+1].next_arg, function_list+last_pos+1);
                function_list[last_pos+1].next_arg = NULL;
            } else if (latex[i+1] == 'y') {
                function_list[last_pos] = new_function(func_extract_y, function_list[last_pos+1].next_arg, function_list+last_pos+1);
                function_list[last_pos+1].next_arg = NULL;
            } else {
                printf("ERROR: invalid character after dot operator\n");
                exit(EXIT_FAILURE);
            }
            i ++;
        }
        // Check if the previous term has a factorial
        else if ((func_pos > 0) && (latex[i] == '!')) {
            shift_blocks(function_list, last_pos, func_pos-last_pos);
            func_pos++;
            function_list[last_pos] = new_function(func_factorial, function_list[last_pos+1].next_arg, function_list+last_pos+1);
            function_list[last_pos+1].next_arg = NULL;
        }
        else if (latex[i] == '\\') {
            cmd_end = i+1;
            while (('a' <= latex[cmd_end]) && (latex[cmd_end] <= 'z')) cmd_end++;
            cmd_len = cmd_end - i - 1;
            cmd_start = i+1;
            arg1_end = 0; arg2_end = 0;
            printf("Found command %.*s, %d\n", cmd_len, latex+cmd_start, cmd_len);
            if ((cmd_len == 4) && (strncmp(latex+cmd_start, "frac", cmd_len) == 0)) {
                arg1_end = extract_braces(latex, cmd_end);
                arg2_end = extract_braces(latex, arg1_end+1);
                arg1_start = cmd_end + 1;
                arg2_start = arg1_end + 2;
                printf("Found fraction: %d-%d over %d-%d\n", arg1_start, arg1_end, arg2_start, arg2_end);
                i = arg2_end;
                oper = func_div;
            } else if ((cmd_len == 6) && (strncmp(latex+cmd_start, "arctan", cmd_len)==0)) {
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply arctangent %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end+6;
                oper = func_arctan;
            } else if ((cmd_len == 2) && (strncmp(latex+cmd_start, "ln", cmd_len)==0)) {
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply ln %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_log;
            } else if ((cmd_len == 3) && (strncmp(latex+cmd_start, "exp", cmd_len)==0)) {
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply exp %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_exp;
            } else if ((cmd_len == 3) && (strncmp(latex+cmd_start, "cos", cmd_len)==0)) {
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply cosine %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_cosine;
            } else if ((cmd_len == 3) && (strncmp(latex+cmd_start, "sin", cmd_len)==0)) {
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply sine %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_sine;
            } else if ((cmd_len == 3) && (strncmp(latex+cmd_start, "det", cmd_len)==0)) {
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply det %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_det;
            } else if ((cmd_len == 3) && (strncmp(latex+cmd_start, "abs", cmd_len)==0)) {
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply abs %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_abs;
            } else if ((cmd_len == 3) && (strncmp(latex+cmd_start, "max", cmd_len)==0)) {
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply max %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_max;
            } else if (((strncmp(latex+cmd_start, "sum", cmd_len) == 0) && (cmd_len == 3)) || ((strncmp(latex+cmd_start, "prod", cmd_len) == 0) && (cmd_len == 4))) {
                // Sums have three parts: a definition, an end point, and an expression to be summed
                // We assume that the entire remaining expression is the sum
                i = cmd_end;
                subscript = extract_braces(latex, i+1);
                superscript = extract_braces(latex, subscript+2);
                printf("Sum found, %d, %d\n", subscript, superscript);
                arg1_end = get_next_match(latex, i+2, subscript, '=');
                printf("Variable is %.*s\n", arg1_end-i-2, latex+i+2);
                printf("Initial value is %.*s\n", subscript-arg1_end-1, latex+arg1_end+1);
                printf("Final value is %.*s\n", superscript-subscript-3, latex+subscript+3);
                last_pos = func_pos;
                if (strncmp(latex+cmd_start, "sum", cmd_len) == 0) function_list[func_pos] = new_function(func_sum, NULL, function_list+func_pos+1);
                else function_list[func_pos] = new_function(func_prod, NULL, function_list+func_pos+1);
                func_pos++;
                function_list[func_pos] = new_value(NULL, 0x40, function_list+func_pos+1);
                func_pos++;
                int prev = func_pos;
                func_pos += PARSE_LATEX_REC(latex+arg1_end+1, subscript-arg1_end-1, function_list+func_pos);
                function_list[prev].next_arg = function_list + func_pos;
                prev = func_pos;
                func_pos += PARSE_LATEX_REC(latex+subscript+3, superscript-subscript-3, function_list+func_pos);
                function_list[prev].next_arg = function_list + func_pos;
                // Create a new variable for the summation index
                strncpy(stringbuf + *string_size, latex+i+2, arg1_end-i-2);
                variable_list[*var_size] = new_variable(stringbuf + *string_size, 0, VARIABLE_IN_SCOPE, NULL);
                *string_size += arg1_end-1;
                *var_size += 1;
                // Parse the expression
                func_pos += PARSE_LATEX_REC(latex+superscript+1, end-superscript-1, function_list+func_pos);
                variable_list[*var_size-1].flags &= ~VARIABLE_IN_SCOPE;
                function_list[last_pos+1].value = variable_list + (*var_size-1);
                printf("parsing %.*s done, func_pos %d\n", end, latex, func_pos);
                return func_pos;
            } else if ((cmd_len == 12) && (strncmp(latex+cmd_start, "operatorname", cmd_len) == 0)) {
                arg1_end = extract_braces(latex, cmd_start+12);
                oper = NULL;
                for (int j=0; j < N_FUNCTIONS; j++) {
                    if ((arg1_end-cmd_start-13 == strlen(function_names[j])) && (strncmp(function_names[j], latex+cmd_start+13, arg1_end-cmd_start-13) == 0)) {
                        oper = function_pointers[j];
                        break;
                    }
                }
                if (!oper) {
                    printf("ERROR: unrecognized function %.*s\n", arg1_end-cmd_start-13, latex+cmd_start+13);
                    exit(EXIT_FAILURE);
                }
                i = arg1_end+1;
                arg1_start = i+6;
                arg1_end = extract_parenthetical(latex, i);
                i = arg1_end;
                arg1_end -= 6;
/*            } else if ((cmd_len == 12) && (cmd_start+17 < end) && (strncmp(latex+cmd_start, "operatorname{mod}", cmd_len+5) == 0)) {
                cmd_end += 5;
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply mod %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_mod;
            } else if ((cmd_len == 12) && (cmd_start+19 < end) && (strncmp(latex+cmd_start, "operatorname{floor}", cmd_len+7) == 0)) {
                cmd_end += 7;
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply floor %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_floor;
            } else if ((cmd_len == 12) && (cmd_start+21 < end) && (strncmp(latex+cmd_start, "operatorname{polygon}", cmd_len+9) == 0)) {
                cmd_end += 9;
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply polygon %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_polygon;
            } else if ((cmd_len == 12) && (cmd_start+19 < end) && (strncmp(latex+cmd_start, "operatorname{total}", cmd_len+7) == 0)) {
                cmd_end += 7;
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply total %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_total;
            } else if ((cmd_len == 12) && (cmd_start+22 < end) && (strncmp(latex+cmd_start, "operatorname{distance}", cmd_len+10) == 0)) {
                cmd_end += 10;
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply floor %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_distance;
            } else if ((cmd_len == 12) && (cmd_start+17 < end) && (strncmp(latex+cmd_start, "operatorname{rgb}", cmd_len+5) == 0)) {
                cmd_end += 5;
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply rgb %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_rgb;
            } else if ((cmd_len == 12) && (cmd_start+17 < end) && (strncmp(latex+cmd_start, "operatorname{max}", cmd_len+5) == 0)) {
                cmd_end += 5;
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply max %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_max;
            } else if ((cmd_len == 12) && (cmd_start+17 < end) && (strncmp(latex+cmd_start, "operatorname{abs}", cmd_len+5) == 0)) {
                cmd_end += 5;
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply abs %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_abs;
            } else if ((cmd_len == 12) && (cmd_start+18 < end) && (strncmp(latex+cmd_start, "operatorname{sort}", cmd_len+6) == 0)) {
                cmd_end += 6;
                arg1_end = extract_parenthetical(latex, cmd_end);
                printf("multiply max %.*s\n", arg1_end - cmd_end - 12, latex+cmd_end+6);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end + 6;
                oper = func_sort;*/
            } else if (((i == 0) || ((i >= 5) && (strncmp(latex+i-5, "\\cdot", 5) == 0))) && (cmd_len == 4) && (strncmp(latex+cmd_start, "left[", cmd_len+1) == 0)) {
                // Brackets can only be interpreted as a list if they are the first thing in the expression
                // or if they follow a \cdot.
                arg1_end = extract_brackets(latex, i);
                n_terms = 0;
                for (int j=i+6; j < end; j++) {
                    if (latex[j] == '{') n_braces++;
                    else if (latex[j] == '}') n_braces--;
                    else if (strncmp(latex+j, "\\left", 5) == 0) {
                        n_leftright++;
                        j += 4;
                    } else if (strncmp(latex+j, "\\right", 6) == 0) {
                        n_leftright--;
                        j += 5;
                    } else if ((j+17 < end) && (strncmp(latex+j, "\\operatorname{for}", 18) == 0) && (n_leftright == 0) && (n_braces == 0)) {
                        // Decompose for expression
                        // These expressions are of the form expr0\operatorname{for}var1=expr1,var2=expr2,:::,varn=exprn
                        // To parse a for expression, one must first create the variable declarations, then parse expr0,
                        // then parse expr1, expr2, and so on.
                        int old_var_size = *var_size;
                        
                        printf("for expression detected: %.*s\n", j-i-6, latex+i+6);
                        // Insert the func_for block
                        last_pos = func_pos;
                        function_list[func_pos] = new_function(func_for, NULL, function_list+func_pos+1);
                        func_pos++;
                        // Parse the list of variable definitions
                        start = j;
                        j += 18;
                        while (1>0) {
                            subscript = get_next_match(latex, j, arg1_end+1, '=');
                            printf("for expression has variable %.*s", subscript-j, latex+j);
                            // Insert the variable into the variable table. pointer will point to the beginning
                            // of the definition of that variable in latex. type will contain the length of
                            // the definition. These will be cleared once the variable definitions have been parsed.
                            // Doing this now ensures the variables follow one another in the variable table
                            strncpy(stringbuf + *string_size, latex+j, subscript-j);
                            variable_list[*var_size] = new_variable(stringbuf + *string_size, 0, VARIABLE_IN_SCOPE, (void*)(latex+subscript+1));
                            *var_size += 1;
                            *string_size += subscript-j+1;
                            j = get_next_match(latex, subscript, arg1_end+1, ',');
                            if (j < 0) {
                                variable_list[*var_size - 1].type = arg1_end-6-subscript-1;
                                printf(" with value expression %p %.*s\n", latex+subscript+1, end-subscript-1, latex+subscript+1);
                                break;
                            }
                            variable_list[*var_size - 1].type = j-subscript-1;
                            printf(" with value expression %.*s\n", j-subscript-1, latex+subscript+1);
                            j++;
                        }
                       
                        // Parse the main expression
                        uint32_t last_term = func_pos;
                        func_pos += PARSE_LATEX_REC(latex+i+6, start-i-6, function_list+func_pos);

                        // Parse the definitions of the variables
                        // Clear the variables and take them out of scope
                        for (j=old_var_size; j < *var_size; j++) {
                            printf("recursing on %p %.*s\n", variable_list[j].pointer, variable_list[j].type, (char*)(variable_list[j].pointer));
                            function_list[last_term].next_arg = function_list+func_pos;
                            last_term = func_pos;
                            function_list[func_pos] = new_function(func_equals, NULL, function_list+func_pos+1);
                            func_pos++;
                            function_list[func_pos] = new_value(variable_list+j, 0x40, function_list+func_pos+1);
                            func_pos++;
                            func_pos += PARSE_LATEX_REC((char*)(variable_list[j].pointer), variable_list[j].type, function_list+func_pos);
                            variable_list[j].flags &= ~VARIABLE_IN_SCOPE;
                            variable_list[j].pointer = NULL;
                            variable_list[j].type = 0;
                        }

                        n_terms++;
                        break;
                    }
                }
                i = arg1_end;
                if (n_terms) {
                    arg1_end = 0;
                    arg1_start = 0;
                } else {
                    arg1_end -= 6;
                    arg1_start = cmd_end+1;
                    oper = func_list;
                }
            } else if ((cmd_len == 4) && (strncmp(latex+cmd_start, "left(", cmd_len+1) == 0)) {
                arg1_end = extract_parenthetical(latex, cmd_start-1);
                i = arg1_end;
                func_pos = insert_product_term(function_list, last_pos, func_pos);
                last_pos = func_pos;
                flags = 0;
                func_pos += PARSE_LATEX_REC(latex+cmd_start+5, arg1_end - cmd_start - 11, function_list+func_pos);
                if (flags & PARSE_ACTION) {
                    shift_blocks(function_list, last_pos, func_pos-last_pos);
                    func_pos++;
                    function_list[last_pos] = new_function(func_chain_actions, NULL, function_list+last_pos+1);
                    *result_flags |= PARSE_ACTION;
                }
                else if (flags & PARSE_COMMA) {
                    printf("parsing point\n");
                    shift_blocks(function_list, last_pos, func_pos-last_pos);
                    func_pos++;
                    function_list[last_pos] = new_function(func_point, NULL, function_list+last_pos+1);
                }
                arg1_end = 0;
                arg2_end = 0;
            } else if ((cmd_len == 4) && (strncmp(latex+cmd_start, "left\\{", cmd_len+2) == 0)) {
                arg1_end = extract_braces(latex, cmd_start+5);
                printf("Found conditional from %d to %d: %.*s\n", cmd_start+6, arg1_end-7, arg1_end-cmd_start-13, latex+cmd_start+6);
                last_pos = func_pos;
                function_list[func_pos] = new_function(func_conditional, NULL, function_list+func_pos+1);
                func_pos++;
                i = cmd_start+5;
                uint8_t chain = 0;
                while (1>0) {
                    arg2_start = get_next_match(latex, i+1, arg1_end-7, ':');
                    if (arg2_start < 0) break;
                    printf("colon found at %d, %p\n", arg2_start, function_list+func_pos);
                    if (chain) function_list[last_pos].next_arg = function_list + func_pos;
                    chain = 1;
                    last_pos = func_pos;
                    flags=0;
                    func_pos += PARSE_LATEX_REC(latex+i+1, arg2_start-i-1, function_list+func_pos);
                    *result_flags |= flags;
                    i = get_next_match(latex, arg2_start, arg1_end-7, ',');
                    if (i < 0) break;
                    printf("comma found at %d, %p\n", i, function_list+func_pos);
                    function_list[last_pos].next_arg = function_list + func_pos;
                    last_pos = func_pos;
                    flags = 0;
                    func_pos += PARSE_LATEX_REC(latex+arg2_start+1, i-arg2_start-1, function_list+func_pos);
                    if ((flags & PARSE_ACTION) && (function_list[last_pos].oper != func_chain_actions)) {
                        // If the parsed result is a single action, it must be wrapped
                        shift_blocks(function_list, last_pos, func_pos-last_pos);
                        func_pos++;
                        function_list[last_pos] = new_function(func_chain_actions, NULL, function_list+last_pos+1);
                    }
                    *result_flags |= flags;
                }
                // Failed to find colon, last separator is comma, ends with a default expression
                if (arg2_start < 0) {
                    function_list[last_pos].next_arg = function_list + func_pos;
                    last_pos = func_pos;
                    flags = 0;
                    func_pos += PARSE_LATEX_REC(latex+i+1, arg1_end-i-8, function_list+func_pos);
                    *result_flags |= flags;
                } 
                // Failed to find comma, last separator is colon, ends with no default expression
                else {
                    function_list[last_pos].next_arg = function_list + func_pos;
                    last_pos = func_pos;
                    func_pos += PARSE_LATEX_REC(latex+arg2_start+1, arg1_end-arg2_start-8, function_list+func_pos);
                    *result_flags |= flags;
                }
                printf("result flags are %02x\n", *result_flags);
                    
                i = arg1_end;
                arg1_end = 0;
                arg2_end = 0;
            } else {
                // Otherwise, interpret as a variable (such as \alpha or \beta)
                int varindex = -1;
                uint8_t flags;
                for (int j=0; j < *var_size; j++) {
                    uint8_t flags = variable_list[j].flags;
                    //printf("comparin %.*s with %.*s, lengths %d and %d\n", cmd_len+1, variable_list[j].name, cmd_len+1, latex+cmd_start-1, strlen(variable_list[j].name)
                    if ((flags & VARIABLE_IN_SCOPE) && (strncmp(variable_list[j].name, latex+cmd_start-1, cmd_len+1) == 0) && 
                            (strlen(variable_list[j].name) == cmd_len+1)) {
                        varindex = j;
                        if (flags & VARIABLE_ARGUMENT) break; // If an argument is in scope, prefer it
                    }
                }
                if (varindex == -1) {
                    for (int j=0; j < *var_size; j++) printf("variable %s at %p has flags %02x\n", variable_list[j].name, variable_list+j, variable_list[j].flags);
                    printf("ERROR: variable %.*s not found!\n", cmd_len+1, latex+cmd_start-1);
                    exit(EXIT_FAILURE);
                }
                if (variable_list[varindex].flags & VARIABLE_FUNCTION) {
                    // Function found
                    // We cannot assume that the function has been defined, only declared. The function
                    // definitions must be correctly connected later. For now, we connect to the variable
                    // declaration
                    arg1_end = extract_parenthetical(latex, cmd_end);
                    printf("multiply %.*s at %d with argument %.*s\n", cmd_len+1, latex+cmd_start-1, varindex, arg1_end-cmd_end+1, latex+cmd_end);
                    func_pos = insert_product_term(function_list, last_pos, func_pos);
                    function_list[func_pos] = new_function(func_user_defined, NULL, function_list+func_pos+1);
                    function_list[func_pos].value = variable_list+varindex;
                    last_pos = func_pos;
                    func_pos++;
                    func_pos += PARSE_LATEX_REC(latex+cmd_end+6, arg1_end-cmd_end-12, function_list+func_pos);
                    i = arg1_end;
                } else {
                    // Variable found
                    printf("multiply %.*s at %d, func_pos %d\n", cmd_len+1, latex+cmd_start-1, varindex, func_pos);
                    func_pos = insert_product_term(function_list, last_pos, func_pos);
                    printf("func_pos %d\n", func_pos);
                    function_list[func_pos] = new_value(variable_list+varindex, 0x40, NULL);
                    last_pos = func_pos;
                    func_pos++;
                    i = cmd_end-1;
                }
                arg1_end = 0;
                arg2_end = 0;
            }
            if (arg2_end) {
                // Function decomposition for two arguments
                func_pos = insert_product_term(function_list, last_pos, func_pos);
                function_list[func_pos] = new_function(oper, NULL, function_list+func_pos+1);
                last_pos = func_pos;
                func_pos++;
                func_pos += PARSE_LATEX_REC(latex+arg1_start, arg1_end - arg1_start, function_list+func_pos);
                function_list[last_pos+1].next_arg = function_list + func_pos;
                func_pos += PARSE_LATEX_REC(latex+arg2_start, arg2_end - arg2_start, function_list+func_pos);
            } else if (arg1_end) {
                // Function decomposition for one argument
                func_pos = insert_product_term(function_list, last_pos, func_pos);
                function_list[func_pos] = new_function(oper, NULL, function_list+func_pos+1);
                last_pos = func_pos;
                func_pos++;
                func_pos += PARSE_LATEX_REC(latex+arg1_start, arg1_end - arg1_start, function_list+func_pos);
                // For empty lists
                if (last_pos+1 == func_pos) function_list[last_pos].first_arg = NULL;
            }
        }
        else if ((('a' <= latex[i]) && (latex[i] <= 'z')) || (('A' <= latex[i]) && (latex[i] <= 'Z'))) {
            subscript = 0; superscript = 0; start = i;
            if ((i+1 < end) && (latex[i+1] == '_')) {
                subscript = extract_braces(latex, i+2);
                i = subscript;
            }
            // We must parse superscripts separately here because functions like sin, cos, tan, etc.
            // can be exponentiated directly
            if ((i+1 < end) && (latex[i+1] == '^')) {
                superscript = extract_braces(latex, i+2);
                i = superscript;
            }
            if (subscript == 0) subscript = start;
            if (superscript) {
                func_pos = insert_product_term(function_list, last_pos, func_pos);
                function_list[func_pos] = new_function(func_exponentiate, NULL, function_list+func_pos+1);
                last_pos = func_pos;
                func_pos++;
                printf("multiply %.*s to the power of %.*s\n", subscript-start+1, latex+start, superscript-subscript-3, latex+subscript+3);
                func_pos += PARSE_LATEX_REC(latex+start, subscript-start+1, function_list+func_pos);
                function_list[last_pos+1].next_arg = function_list + func_pos;
                func_pos += PARSE_LATEX_REC(latex+subscript+3, superscript-subscript-3, function_list+func_pos);
            } else {
                // Multiply by a single variable or function
                int varindex = -1;
                uint8_t flags;
                for (int j=0; j < *var_size; j++) {
                    uint8_t flags = variable_list[j].flags;
                    if ((flags & VARIABLE_IN_SCOPE) && (strncmp(variable_list[j].name, latex+start, subscript-start+1) == 0) && (strlen(variable_list[j].name)==subscript-start+1)) {
                        varindex = j;
                        if (flags & VARIABLE_ARGUMENT) break; // If an argument is in scope, prefer it
                    }
                }
                if (varindex == -1) {
                    // If the variable is not found, we may have to create one. The variable will remain
                    // in-scope until parsing is done. The PARSE_NEWVAR flag is set
                    if ((*var_size > 0) && (variable_list[*var_size-1].flags & VARIABLE_IN_SCOPE) && (variable_list[*var_size-1].flags & VARIABLE_XYLIKE)){
                        printf("ERROR: variable %.*s not found!\n", subscript-start+1, latex+start);
                        exit(EXIT_FAILURE);
                    }
                    printf("Creating variable %.*s\n", subscript-start+1, latex+start);
                    strncpy(stringbuf + *string_size, latex+start, subscript-start+1);
                    variable_list[*var_size] = new_variable(stringbuf + *string_size, 0, VARIABLE_XLIKE | VARIABLE_IN_SCOPE, NULL);
                    varindex = *var_size;
                    *var_size += 1;
                    *string_size += subscript-start+2;
                }
                if (variable_list[varindex].flags & VARIABLE_FUNCTION) {
                    // Function found
                    // We cannot assume that the function has been defined, only declared. The function
                    // definitions must be correctly connected later. For now, we connect to the variable
                    // declaration
                    arg1_end = extract_parenthetical(latex, subscript+1);
                    printf("multiply %.*s at %d with argument %.*s\n", subscript-start+1, latex+start, varindex, arg1_end-subscript, latex+subscript+1);
                    func_pos = insert_product_term(function_list, last_pos, func_pos);
                    function_list[func_pos] = new_function(func_user_defined, NULL, function_list+func_pos+1);
                    function_list[func_pos].value = variable_list+varindex;
                    last_pos = func_pos;
                    func_pos++;
                    func_pos += PARSE_LATEX_REC(latex+subscript+7, arg1_end-subscript-13, function_list+func_pos);
                    i = arg1_end;
                } else {
                    // Variable found
                    printf("multiply %.*s at %d, func_pos %d\n", subscript-start+1, latex+start, varindex, func_pos);
                    func_pos = insert_product_term(function_list, last_pos, func_pos);
                    printf("func_pos %d\n", func_pos);
                    function_list[func_pos] = new_value(variable_list+varindex, 0x40, NULL);
                    last_pos = func_pos;
                    func_pos++;
                }
            }
        } else if (('0' <= latex[i]) && (latex[i] <= '9')) {
            i = extract_double(latex, i, &value);
            func_pos = insert_product_term(function_list, last_pos, func_pos);
            printf("multiply constant %f, stack %d, func_pos %d\n", value, *stack_size, func_pos);
            stack[*stack_size] = value;
            function_list[func_pos] = new_value(stack + (*stack_size), 1<<8, NULL);
            last_pos = func_pos;
            func_pos++;
            *stack_size += 1;
        }

        // Check if the term is followed by a \cdot. If this is the case, we skip the index check.
        // This allows the following bracketed term to be interpreted as a list rather than an
        // index to the previous term.
        if ((i+5 < end) && (strncmp(latex+i+1, "\\cdot", 5) == 0)) {
            i += 5;
        }


    }
    
    
    printf("parsing %.*s done, func_pos %d\n", end, latex, func_pos);
    return func_pos;
}

/*
00 0x103ba7200: first_arg: 0x103ba7228, next_arg: 0x0, oper: 0x103aea100                    (add)
01 0x103ba7228: first_arg: 0x103ba7250, next_arg: 0x103ba72a0, oper: 0x103aea200            (mult)
02 0x103ba7250: first_arg: 0x0, next_arg: 0x103ba7278, oper: 0x0                            2
03 0x103ba7278: first_arg: 0x0, next_arg: 0x0, oper: 0x0                                    x
04 0x103ba72a0: first_arg: 0x103ba72c8, next_arg: 0x103ba7390, oper: 0x103aea200            (mult)
05 0x103ba72c8: first_arg: 0x103ba72f0, next_arg: 0x0, oper: 0x103aea300                    (exp)
06 0x103ba72f0: first_arg: 0x0, next_arg: 0x103ba7318, oper: 0x0                            x
07 0x103ba7318: first_arg: 0x103ba7340, next_arg: 0x0, oper: 0x103aea100                    (add)
08 0x103ba7340: first_arg: 0x0, next_arg: 0x103ba7368, oper: 0x0                            4
09 0x103ba7368: first_arg: 0x0, next_arg: 0x0, oper: 0x0                                    x
10 0x103ba7390: first_arg: 0x103ba73b8, next_arg: 0x4008000000000000, oper: 0x103aea100     add
11 0x103ba73b8: first_arg: 0x0, next_arg: 0x103ba73e0, oper: 0x0                            7
12 0x103ba73e0: first_arg: 0x0, next_arg: 0x0, oper: 0x0                                    3


00 0x10d120200: first_arg: 0x10d120228, next_arg: 0x0, oper: 0x10d062f70                            (add)
01 0x10d120228: first_arg: 0x10d120250, next_arg: 0x10d1202a0, oper: 0x10d063070                    (mult)
02 0x10d120250: first_arg: 0x0, next_arg: 0x10d120278, oper: 0x0                                    2
03 0x10d120278: first_arg: 0x0, next_arg: 0x0, oper: 0x0                                            x
04 0x10d1202a0: first_arg: 0x10d1202c8, next_arg: 0x10d120430, oper: 0x10d063070                    (mult)
05 0x10d1202c8: first_arg: 0x10d1202f0, next_arg: 0x10d120390, oper: 0x10d063170                    (exp)
06 0x10d1202f0: first_arg: 0x0, next_arg: 0x10d120318, oper: 0x0                                    x
07 0x10d120318: first_arg: 0x10d120340, next_arg: 0x0, oper: 0x10d062f70                            (add)
08 0x10d120340: first_arg: 0x0, next_arg: 0x10d120368, oper: 0x0                                    4
09 0x10d120368: first_arg: 0x0, next_arg: 0x0, oper: 0x0                                            x
10 0x10d120390: first_arg: 0x4008000000000000, next_arg: 0x401c000000000000, oper: 0x10d062d50      (cos)
11 0x10d1203b8: first_arg: 0x10d1203e0, next_arg: 0x0, oper: 0x10d062f70                            (add)
12 0x10d1203e0: first_arg: 0x0, next_arg: 0x10d120408, oper: 0x0                                    4
13 0x10d120408: first_arg: 0x0, next_arg: 0x0, oper: 0x0                                            x
14 0x10d120430: first_arg: 0x10d120458, next_arg: 0x0, oper: 0x10d062f70                            (add)
15 0x10d120458: first_arg: 0x0, next_arg: 0x10d120480, oper: 0x0                                    7
16 0x10d120480: first_arg: 0x0, next_arg: 0x0, oper: 0x0                                            3
*/

void parse_latex(char *latex, function *function_list, double *stack, variable *variable_list, char *stringbuf, int *stack_size, int *var_size, int *string_size) {
    // Find all of the variables in the expression
    int end = strlen(latex);
    int i=0;
    // Parse LaTeX
    uint8_t flags = 0;
    int func_pos = parse_latex_rec(latex, end, function_list, stack, variable_list, stringbuf, stack_size, var_size, string_size, &flags);

    for (i=0; i < func_pos; i++) {
        printf("%02d %p: first_arg: %p, next_arg: %p, oper: %p, value: %p\n", i, function_list+i, function_list[i].first_arg, function_list[i].next_arg, function_list[i].oper, function_list[i].value);
    }

    for (i=0; i < *var_size; i++)
        printf("%s variable pointer at %p, points to %p, type %08X\n", variable_list[i].name, &(variable_list[i].pointer), variable_list[i].pointer, variable_list[i].type);
}

expression* topological_sort_rec(expression *start) {
    expression *top = start;
    // Iterate over the dependencies and decrease the number of dependents
    for (int i=0; i < start->num_dependencies; i++) {
        (start->dependencies)[i]->num_dependents--;
        if ((start->dependencies)[i]->num_dependents == 0) {
            // If a dependency is found that no longer depends on anything
            // that has not been added to the chain, recurse on that
            // dependency
            (start->dependencies)[i]->next_expr = top;
            top = topological_sort_rec((start->dependencies)[i]);
        }
    }
    return top;
}

expression* topological_sort(expression *expression_list, int n_expr) {
    expression *top = NULL;
    for (int i=0; i < n_expr; i++) {
        if ((expression_list[i].num_dependents == 0) && (expression_list[i].next_expr == NULL)) {
            // Node may be used as a starting point
            expression_list[i].next_expr = top;
            top = topological_sort_rec(expression_list+i);
        }
    }
    return top;
}

uint8_t varncmp(char *varname, char *target, uint32_t len) {
    return (strncmp(varname, target, len) == 0) && (strlen(varname) == len);
}

void insert_interval_functions(function *func) {
    const oper_data *oper = oper_lookup(func->oper);
    func->inter = oper->inter;
    if (!(func->inter)) return;
    function *arg = func->first_arg;
    while (arg) {
        insert_interval_functions(arg);
        if (!(arg->inter)) {
            func->inter = NULL;
            return;
        }
        arg = arg->next_arg;
    }
}

void load_file(char *fname, file_data *fd) {
    expression *expression_list = fd->expression_list;
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
        printf("ERROR: unable to open file %s\n", fname);
        exit(EXIT_FAILURE);
    }
    uint32_t n_expr = 0;
    char *temp;
    char *line=NULL;
    ssize_t read;
    size_t len=0;
    while ((read = getline(&line, &len, fp)) != -1) {
        if (read <= 1) continue;
        if (line[0] == '#') continue;
        if (line[read-1] == '\n') {
            line[read-1] = 0;
            read--;
        }
        printf("Loaded expression %s\n", line);
        temp = malloc((read+1)*sizeof(char));
        strcpy(temp, line);
        expression_list[n_expr] = new_expression(NULL, NULL, NULL, 0, 0);
        expression_list[n_expr].def = temp;
        if ((n_expr%7)&0x01) expression_list[n_expr].color[0] = 255;
        else expression_list[n_expr].color[0] = 0;
        if ((n_expr%7)&0x02) expression_list[n_expr].color[1] = 255;
        else expression_list[n_expr].color[1] = 0;
        if ((n_expr%7)&0x04) expression_list[n_expr].color[2] = 255;
        else expression_list[n_expr].color[2] = 0;
        n_expr++;
    }
    fd->n_expr = n_expr;
}

//expression* parse_file(function *function_list, double *stack, variable *variable_list, char *stringbuf, expression *expression_list, uint32_t *n_func, uint32_t *n_var, uint32_t n_expr, uint32_t *n_stack) {
expression *parse_file(file_data *fd, char *stringbuf) {
    expression *expression_list = fd->expression_list;
    variable *variable_list = fd->variable_list;
    function *function_list = fd->function_list;
    double *stack = fd->stack;
    uint32_t n_expr = fd->n_expr, *n_var = &(fd->n_var), *n_func = &(fd->n_func), *n_stack = &(fd->n_stack);
    size_t len=0;
    ssize_t read;
    int i;
    int stack_size = 0;
    variable *varpos = variable_list;
    expression *exprpos = expression_list;
    while (varpos[0].name) varpos++;
    char *stringpos = stringbuf;
    int varidx;
    uint8_t in_command;
    uint8_t var_alpha;
    uint8_t var_beta;
    uint8_t flags=0;
    printf("variable_list: %p\nexpression_list: %p\nstack: %p\nstringpos: %p\nn_expr: %d\n", variable_list, expression_list, stack, stringpos, n_expr);

    char *line;
    for (uint32_t expr_idx = 0; expr_idx < n_expr; expr_idx++) {
        line = expression_list[expr_idx].def;
        read = strlen(line);
        exprpos = expression_list+expr_idx;

        printf("Found expression of length %zd\n", read);
        printf("    %s\n", line);
        i = 1;
        var_alpha = ((read >= 6) && (strncmp(line, "\\alpha", 6) == 0));
        var_beta = ((read >= 5) && (strncmp(line, "\\beta", 5) == 0));
        if (var_alpha) i += 5;
        if (var_beta) i += 4;
        if ((('a' <= line[0]) && (line[0] <= 'z')) || (('A' <= line[0]) && (line[0] <= 'Z')) || var_alpha || var_beta) {
            // If the line starts with a letter, it could be a variable or function definition.
            // First, we need to look for a subscript and determine the range in which the
            // variable or function name lies.
            if (line[i] == '_') {
                i = extract_braces(line, i+1)+1;
            }
            // Check if the line could be a function definition
            if ((i+6 < read) && (strncmp(line+i, "\\left(", 6) == 0)) {
                int arg = extract_parenthetical(line, i);
                // If no equals sign is found, then it is not a definition. In that case, continue
                // and increment the expression list counter
                if (line[arg+1] == '=') printf("    Found function %.*s with args %.*s\n", i, line, arg+1-i, line+i);
                else continue;
                
                // Check if the function has already been defined
                for (int j=0; variable_list+j < varpos; j++) {
                    if (strcmp(variable_list[j].name, line) == 0) {
                        printf("ERROR: multiple definition for %.*s\n", i, line);
                        exit(EXIT_FAILURE);
                    }
                }
                strncpy(stringpos, line, i);
                // Pointer is NULL because the function has not been defined yet
                varpos[0] = new_variable(stringpos, 0, VARIABLE_IN_SCOPE | VARIABLE_FUNCTION, NULL);
                exprpos[0].var = varpos;
                stringpos += i+1;
                varpos++;

                // Parse the arguments of the function and add the arguments of the function
                // to the variable list.
                for (int j=i+6; j < arg-6; j++) {
                    if (line[j] == ',') continue;
                    if (line[j] == ' ') continue;
                    if ((j < arg-7) && (line[j+1] == '_')) i = extract_braces(line, j+2);
                    else i = j;
                    strncpy(stringpos, line+j, i+1-j);
                    printf("        argument %.*s\n", i+1-j, line+j);
                    varpos[0] = new_variable(stringpos, 0, VARIABLE_ARGUMENT, NULL);
                    stringpos += i+2-j;
                    varpos++;
                    j=i;
                }
                i = arg+2;
            } else if (line[i] == '=') {
                // If an equals sign immediately follows the variable name, then it is a variable definition
                printf("    Found variable %.*s\n", i, line);
                //int end = extract_double(line, i+1, stack+stack_size);
                strncpy(stringpos, line, i);
                //if (end+2 == read) {
                //    printf("    Variable has value %f\n", stack[stack_size]);
                //    varpos[0] = new_variable(stringpos, 1<<8, VARIABLE_IN_SCOPE, stack+stack_size);
                //    stack_size++;
                //} else {
                    varpos[0] = new_variable(stringpos, 0, VARIABLE_IN_SCOPE, NULL);
                //}
                exprpos->var = varpos;
                varpos++;
                stringpos += i+1;
                i++;
            } else continue;
            // In either case (variable or function definition, we need to finish up the expression struct.
            // This involves extracting the variables used in the definition. However, we can't do this
            // until we parse all of the expressions. Thus, we store the index to the expression struct.
            exprpos->expr_begin = i;

        }
    }
    printf("loop done, %p, %p\n", variable_list->name, stringpos);
    for (i=0; variable_list+i < varpos; i++)
        printf("%s variable pointer at %p, points to %p, type %08X\n", variable_list[i].name, variable_list+i, variable_list[i].pointer, variable_list[i].type);

    printf("A total of %d expressions found\n", n_expr);
    exprpos = expression_list;
    int subscript;
    expression *deptable[300];
    int deptable_ofs = 0;
    uint8_t already_exists;
    // Store the beginning of the implicit variable section
    variable *old_varpos = varpos;
    int func_pos=0, last_pos=0;
    // Pointer to any new variable created as the input to a particular expression.
    // For example, u in w=10^{u}
    variable *local_variable = NULL;
    for (uint32_t expr_idx=0; expr_idx < n_expr; expr_idx++) {
        line = expression_list[expr_idx].def;
        read = strlen(line);
        exprpos = expression_list + expr_idx;
        
        i = exprpos->expr_begin;
        exprpos->dependencies = deptable+deptable_ofs;
        exprpos->flags |= EXPRESSION_FIXED;
        in_command = 0;
        printf("Checking expression %p\n", exprpos);
        local_variable = NULL;

        for (int j=i; j < read; j++) {
            if (line[j] == '\\') in_command = 1;
            else if ((('a' <= line[j]) && (line[j] <= 'z')) || (('A' <= line[j]) && (line[j] <= 'Z'))) {
                var_alpha = ((j+4 < read) && (strncmp(line+j, "alpha", 5) == 0));
                var_beta = ((j+3 < read) && (strncmp(line+j, "beta", 4)==0));
                if (!in_command || var_alpha || var_beta) {
                    subscript = j;
                    if (var_alpha) subscript += 4;
                    if (var_beta) subscript += 3;
                    if (var_alpha || var_beta) j--;
                    if ((subscript+1 < read) && (line[subscript+1] == '_')) {
                        subscript = extract_braces(line, subscript+2)+1;
                    } else {
                        subscript = subscript+1;
                    }
                    printf("   Found variable %.*s\n", subscript-j, line+j);
                    // Now, we must look throught the expression list and find each variable so we can add it
                    // to the dependency table. First, we check if it hasn't already been added to the table
                    already_exists = 0;
                    for (int k=0; k < exprpos->num_dependencies; k++) {
                        if (varncmp((((exprpos->dependencies)[k])->var)->name, line+j, subscript-j)) {
                            printf("    Variable already in dependency table (%s, %d, %.*s)\n", exprpos->dependencies[k]->var->name, subscript-j, subscript-j, line+j);
                            already_exists = 1;
                            break;
                        }
                    }
                    if (!already_exists) {
                        // If the dependency has not already been added, search for the variable in the
                        // expression table. Note that this will skip variables that are not defined by the
                        // user. For example, it will miss x and y as well as function arguments
                        for (int k=0; k < n_expr; k++) {
                            if ((expression_list[k].var) && (varncmp((expression_list[k].var)->name, line+j, subscript-j))) {
                                printf("    found at position %d\n", k);
                                deptable[deptable_ofs] = expression_list+k;
                                deptable_ofs++;
                                exprpos->num_dependencies++;
                                already_exists = 1;
                                break;
                            }
                        }
                        if (!already_exists && (exprpos->var) && (exprpos->var->flags & VARIABLE_FUNCTION)) {
                            // Iterate over the arguments of the function to check whether the variable is
                            // any of those
                            variable *arg = (exprpos->var)+1;
                            while (arg->flags & VARIABLE_ARGUMENT) {
                                if (varncmp(arg->name, line+j, subscript-j)) {
                                    already_exists = 1;
                                    break;
                                }
                                arg++;
                            }
                        }
                        if (!already_exists && ((line[j] == 'x') || (line[j] == 'y'))) {
                            exprpos->flags |= EXPRESSION_PLOTTABLE;
                            exprpos->flags &= ~EXPRESSION_FIXED;
                        }
                    }
                    
                    j = subscript-1;
                }
            } else if (line[j] == '.') in_command = 1;
            else if (line[j] == '{') {
                if ((j >= 13) && (strncmp(line+j-13, "\\operatorname", 13) == 0)) in_command = 1;
                else in_command = 0;
            } else in_command = 0;
        }
        
        if ((exprpos->var) && (exprpos->var->flags & VARIABLE_FUNCTION)) {
            exprpos->func = function_list+func_pos;
            exprpos->var->pointer = (double*)(function_list+func_pos);
            // Bring the arguments of the function into the local scope for parsing
            variable *first_arg = (exprpos->var)+1;
            while (first_arg->flags & VARIABLE_ARGUMENT) {
                first_arg->flags |= VARIABLE_IN_SCOPE;
                first_arg++;
            }
            // Parse the definition of the function
            last_pos = func_pos;
            int var_size = varpos - variable_list;
            int string_size = stringpos - stringbuf;
            func_pos += parse_latex_rec(line+i, read-i, function_list+func_pos, stack, variable_list, stringbuf, &stack_size, &var_size, &string_size, &flags);
            for (int p=last_pos; p < func_pos; p++) {
                if (function_list[p].oper == func_assign) {
                    exprpos->flags |= EXPRESSION_ACTION;
                    break;
                }
            }
            if (exprpos->flags & EXPRESSION_ACTION) {
                shift_blocks(function_list, last_pos, func_pos-last_pos);
                func_pos++;
                function_list[last_pos] = new_function(func_chain_actions, NULL, function_list+last_pos+1);
                exprpos->var->flags |= VARIABLE_ACTION;
            }

            varpos = var_size + variable_list;
            stringpos = stringbuf + string_size;
            function_list[last_pos].value_type |= TYPE_ABSOLUTE_ADDR;
            // Connect next_arg of the first block in the definition to the variable block of the first argument.
            function_list[last_pos].next_arg = (function*)((exprpos->var)+1);
            // Take the arguments out of the local scope after parsing
            first_arg = (exprpos->var)+1;
            while (first_arg->flags & VARIABLE_ARGUMENT) {
                first_arg->flags &= ~VARIABLE_IN_SCOPE;
                first_arg++;
            }
        } else if (!(exprpos->var) || (exprpos->var->pointer == 0)) {
            exprpos->func = function_list+func_pos;
            int var_size = varpos - variable_list;
            int string_size = stringpos - stringbuf;
            last_pos = func_pos;
            func_pos += parse_latex_rec(line+i, read-i, function_list+func_pos, stack, variable_list, stringbuf, &stack_size, &var_size, &string_size, &flags);
            for (int p=last_pos; p < func_pos; p++) {
                if (function_list[p].oper == func_assign) {
                    exprpos->flags |= EXPRESSION_ACTION;
                    break;
                }
            }
            if (exprpos->flags & EXPRESSION_ACTION) {
                shift_blocks(function_list, last_pos, func_pos-last_pos);
                func_pos++;
                function_list[last_pos] = new_function(func_chain_actions, NULL, function_list+last_pos+1);
                if (exprpos->var) {
                    exprpos->var->flags |= VARIABLE_ACTION;
                    exprpos->var->pointer = (double*)(function_list+last_pos);
                }
            }
            uint8_t no_variables = 1;
            for (int tpos=last_pos; tpos < func_pos; tpos++) {
                if (((function_list[tpos].value_type & 0x40) && (function_list[tpos].oper == func_value)) || (function_list[tpos].oper == func_user_defined)) {
                    no_variables = 0;
                    break;
                }
            }
            if (no_variables) {
                printf("Expression %p has no variables, evaluating\n", exprpos);
                uint32_t argtype = function_list[last_pos].oper(function_list+last_pos, stack+stack_size);
                double *temp = malloc((argtype>>8)*sizeof(double));
                memcpy(temp, stack+stack_size, (argtype>>8)*sizeof(double));
                function_list[last_pos] = new_value(temp, argtype, NULL);
                memset(function_list+last_pos+1, 0, (func_pos-last_pos-1)*sizeof(function));
                func_pos = last_pos+1;
            }
            varpos = var_size + variable_list;
            stringpos = stringbuf + string_size;
            if ((var_size > 0) && (variable_list[var_size-1].flags & VARIABLE_XYLIKE) && (variable_list[var_size-1].flags & VARIABLE_IN_SCOPE)) {
                // If the parser has created a new x or y-like variable, classify the expression as plottable and not fixed
                variable_list[var_size-1].flags &= ~VARIABLE_IN_SCOPE;
                exprpos->flags &= ~EXPRESSION_FIXED;
                exprpos->flags |= EXPRESSION_PLOTTABLE;
            }
        }
        
        if (local_variable) local_variable->flags &= (~VARIABLE_IN_SCOPE);
        
        free(exprpos->def);
    }

    // Connect all of the function calls to their definitions instead of their declarations
    // Connect implicit variables to their actual value (for example, replace u in w=10^{u} with x)
    for (i=0; i < func_pos; i++) {
        if (function_list[i].oper == func_user_defined) {
            function_list[i].value = ((variable*)(function_list[i].value))->pointer;
        }
        if (function_list[i].value_type & 0x40) {
            if (((variable*)(function_list[i].value))->flags & VARIABLE_XLIKE) {
                printf("function block %p used to point to %s (%p), now points to x (%p)\n", function_list+i, ((variable*)(function_list[i].value))->name, function_list[i].value, variable_list);
                function_list[i].value = (void*)variable_list;
            }
            else if (((variable*)(function_list[i].value))->flags & VARIABLE_YLIKE) function_list[i].value = (void*)(variable_list+1);
        }
    }

    // Print the expression table
    for (i=0; i < n_expr; i++) {
        printf("expression pointer at %p, function %p, variable %p, flags %02x, begin %d, dep %p, num %d\n", expression_list+i, expression_list[i].func, expression_list[i].var, expression_list[i].flags, expression_list[i].expr_begin, expression_list[i].dependencies, expression_list[i].num_dependencies);
        insert_interval_functions(expression_list[i].func);
    }
    
    // Print the function table
    for (i=0; i < func_pos; i++) {
        if (function_list[i].value_type & 0x40) printf("%02d %p: first_arg: %p, next_arg: %p, oper: %p, value: %p (%s), value_type: %08x\n", i, function_list+i,
                function_list[i].first_arg, function_list[i].next_arg, function_list[i].oper, function_list[i].value, ((variable*)(function_list[i].value))->name,
                function_list[i].value_type);
        else if (function_list[i].oper == func_value) printf("%02d %p: first_arg: %p, next_arg: %p, oper: %p, value: %p (%f), value_type: %08x\n", i, function_list+i,
                function_list[i].first_arg, function_list[i].next_arg, function_list[i].oper, function_list[i].value, ((double*)(function_list[i].value))[0],
                function_list[i].value_type);
        else printf("%02d %p: first_arg: %p, next_arg: %p, oper: %p, value: %p, value_type: %08x\n", i, function_list+i, function_list[i].first_arg, 
                function_list[i].next_arg, function_list[i].oper, function_list[i].value, function_list[i].value_type);
    }

    FILE *fp = fopen("/tmp/function_list", "w");
    fwrite(&function_list, sizeof(int64_t), 1, fp);
    fwrite(function_list, sizeof(function), func_pos, fp);
    fclose(fp);
    

    // Use topological sort to determine the order in which variable and function assignments must be evaluated.
    for (i=0; i < deptable_ofs; i++) 
        deptable[i]->num_dependents++;
    for (i=0; i < n_expr; i++) {
        if (expression_list[i].num_dependencies > 0) {
            printf("%p (%03d)\n", expression_list+i, i+1);
            for (int j=0; j < expression_list[i].num_dependencies; j++) printf("    %p (%03ld)\n", expression_list[i].dependencies[j], expression_list[i].dependencies[j]-expression_list+1);
        }
    }
    for (i=0; i < n_expr; i++) {
        printf("%d expressions depend on %p\n", expression_list[i].num_dependents, expression_list+i);
    }
    expression *top_expr = topological_sort(expression_list, n_expr);
    printf("Top is %p\n", top_expr);
    expression *expr = top_expr;
    i=0;
    while (expr && (i < n_expr+1)) {
        printf("'%03ld', ", expr-expression_list+1);
        expr = expr->next_expr;
        i++;
    }
    if (i >= n_expr+1) {
        printf("Circular reference!\n");
        exit(EXIT_FAILURE);
    }
    printf("\n");
    
    // Evaluate all variables in accordance with the topological ordering
    printf("stack_size is %d\n", stack_size);
    for (int p=0; p < stack_size; p++) {
        printf("%f\t", stack[p]);
        if (p % 8 == 7) printf("\n");
    }
    printf("\n");
    fd->n_stack = stack_size;
    evaluate_from(fd, top_expr);


    for (i=0; variable_list+i < varpos; i++)
        printf("%s variable pointer at %p, points to %p, type %08X, flags %02X, new_pointer %p\n", variable_list[i].name, variable_list+i, variable_list[i].pointer, variable_list[i].type, variable_list[i].flags, variable_list[i].new_pointer);

    *n_func = func_pos;
    *n_var = varpos - variable_list;

    return top_expr;
}

