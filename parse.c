#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "parse.h"
#include "functions.h"

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
    block.oper = NULL;
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

int parse_latex_rec(char *latex, int end, function *function_list, double *stack, variable *variable_list, int *stack_size, uint8_t *result_flags) {
    printf("Parsing %.*s\n", end, latex);
    
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

    // Decompose commas
    int last_term = 0;
    uint8_t subtract = 0;
    uint8_t start_paren = 0;
    uint8_t n_leftright = 0;
    uint8_t n_braces = 0;
    uint8_t n_terms = 0;
    int func_pos = 0;
    int last_pos = 0;
    uint8_t flags;
    for (int i=0; i < end; i++) {
        if (latex[i] == '{') n_braces++;
        else if (latex[i] == '}') n_braces--;
        else if (strncmp(latex+i, "\\left", 5) == 0) {
            n_leftright++;
            if (((i == 0) || (latex[i-1] == ',')) && (latex[i+5] == '(')) start_paren = 1;
            i += 4;
        } else if (strncmp(latex+i, "\\right", 6) == 0) {
            n_leftright--;
            if (!((i+7 >= end) || (latex[i+7] == ',')) && (latex[i+6] == ')') && (n_leftright == 0)) start_paren = 0;
            i += 5;
        } else if ((latex[i] == ',') && (n_leftright == 0) && (n_braces == 0)) {
            printf("comma term %.*s\n", i - last_term, latex+last_term);
            
            if (n_terms > 0) function_list[last_pos].next_arg = function_list + func_pos;
            last_pos = func_pos;
            
            if (start_paren) {
                func_pos += parse_latex_rec(latex+last_term+6, i - last_term - 13, function_list+func_pos, stack, variable_list, stack_size, &flags);
            } else {
                func_pos += parse_latex_rec(latex+last_term, i - last_term, function_list+func_pos, stack, variable_list, stack_size, &flags);
            }
            last_term = i+1;
            n_terms++;
        }
    }
    if (n_terms) {
        *result_flags |= PARSE_COMMA;
        printf("comma term (final) %.*s\n", end - last_term, latex+last_term);
        
        function_list[last_pos].next_arg = function_list + func_pos;
        last_pos = func_pos;
            
        if (start_paren) {
            func_pos += parse_latex_rec(latex+last_term+6, end - last_term - 13, function_list+func_pos, stack, variable_list, stack_size, &flags);
        } else {
            func_pos += parse_latex_rec(latex+last_term, end - last_term, function_list+func_pos, stack, variable_list, stack_size, &flags);
        }
        if (subtract) function_list[last_pos].value_type ^= 0x80;
        printf("parsing %.*s done, func_pos %d\n", end, latex, func_pos);
        return func_pos;
    }
    
    // Decompose addition
    n_terms = 0;
    for (int i=0; i < end; i++) {
        if (latex[i] == '{') n_braces++;
        else if (latex[i] == '}') n_braces--;
        else if (strncmp(latex+i, "\\left", 5) == 0) {
            n_leftright++;
            if (((i == 0) || (latex[i-1] == '-') || (latex[i-1] == '+')) && (latex[i+5] == '(')) start_paren = 1;
            i += 4;
        } else if (strncmp(latex+i, "\\right", 6) == 0) {
            n_leftright--;
            if (!((i+7 >= end) || (latex[i+7] == '-') || (latex[i+7] == '+')) && (latex[i+6] == ')') && (n_leftright == 0)) start_paren = 0;
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
            
            if (start_paren) {
                func_pos += parse_latex_rec(latex+last_term+6, i - last_term - 13, function_list+func_pos, stack, variable_list, stack_size, &flags);
            } else {
                func_pos += parse_latex_rec(latex+last_term, i - last_term, function_list+func_pos, stack, variable_list, stack_size, &flags);
            }
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
            
        if (start_paren) {
            func_pos += parse_latex_rec(latex+last_term+6, end - last_term - 13, function_list+func_pos, stack, variable_list, stack_size, &flags);
        } else {
            func_pos += parse_latex_rec(latex+last_term, end - last_term, function_list+func_pos, stack, variable_list, stack_size, &flags);
        }
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
    int subscript, superscript, start;
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
            uint32_t parsed_length = parse_latex_rec(latex+i+6, idx_end - i - 12, function_list+func_pos, stack, variable_list, stack_size, &flags);
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
            func_pos += parse_latex_rec(latex+2+i, superscript - i - 2, function_list+func_pos, stack, variable_list, stack_size, &flags);
            i = superscript;
        } else if (latex[i] == '\\') {
            cmd_end = i+1;
            while (('a' <= latex[cmd_end]) && (latex[cmd_end] <= 'z')) cmd_end++;
            cmd_len = cmd_end - i - 1;
            cmd_start = i+1;
            arg1_end = 0; arg2_end = 0;
            printf("Found command %.*s\n", cmd_len, latex+cmd_start);
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
            } else if (((i == 0) || ((i >= 5) && (strncmp(latex+i-5, "\\cdot", 5) == 0))) && (cmd_len == 4) && (strncmp(latex+cmd_start, "left[", cmd_len+1) == 0)) {
                // Brackets can only be interpreted as a list if they are the first thing in the expression
                // or if they follow a \cdot.
                arg1_end = extract_brackets(latex, i);
                i = arg1_end;
                arg1_end -= 6;
                arg1_start = cmd_end+1;
                oper = func_list;
            }
            if ((cmd_len == 4) && (strncmp(latex+cmd_start, "left(", cmd_len+1) == 0)) {
                arg1_end = extract_parenthetical(latex, cmd_start-1);
                i = arg1_end;
                func_pos = insert_product_term(function_list, last_pos, func_pos);
                last_pos = func_pos;
                flags = 0;
                func_pos += parse_latex_rec(latex+cmd_start+5, arg1_end - cmd_start - 11, function_list+func_pos, stack, variable_list, stack_size, &flags);
                if (flags & PARSE_COMMA) {
                    printf("parsing point\n");
                    shift_blocks(function_list, last_pos, func_pos-last_pos);
                    func_pos++;
                    function_list[last_pos] = new_function(func_point, NULL, function_list+last_pos+1);
                }
            } else if (arg2_end) {
                // Function decomposition for two arguments
                func_pos = insert_product_term(function_list, last_pos, func_pos);
                function_list[func_pos] = new_function(oper, NULL, function_list+func_pos+1);
                last_pos = func_pos;
                func_pos++;
                func_pos += parse_latex_rec(latex+arg1_start, arg1_end - arg1_start, function_list+func_pos, stack, variable_list, stack_size, &flags);
                function_list[last_pos+1].next_arg = function_list + func_pos;
                func_pos += parse_latex_rec(latex+arg2_start, arg2_end - arg2_start, function_list+func_pos, stack, variable_list, stack_size, &flags);
            } else {
                // Function decomposition for one argument
                func_pos = insert_product_term(function_list, last_pos, func_pos);
                function_list[func_pos] = new_function(oper, NULL, function_list+func_pos+1);
                last_pos = func_pos;
                func_pos++;
                func_pos += parse_latex_rec(latex+arg1_start, arg1_end - arg1_start, function_list+func_pos, stack, variable_list, stack_size, &flags);
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
                func_pos += parse_latex_rec(latex+start, subscript-start+1, function_list+func_pos, stack, variable_list, stack_size, &flags);
                function_list[last_pos+1].next_arg = function_list + func_pos;
                func_pos += parse_latex_rec(latex+subscript+3, superscript-subscript-3, function_list+func_pos, stack, variable_list, stack_size, &flags);
            } else {
                // Multiply by a single variable or function
                int varindex = -1;
                uint8_t flags;
                for (int j=0; variable_list[j].name; j++) {
                    uint8_t flags = variable_list[j].flags;
                    if ((flags & VARIABLE_IN_SCOPE) && (strncmp(variable_list[j].name, latex+start, subscript-start+1) == 0)) {
                        varindex = j;
                        if (flags & VARIABLE_ARGUMENT) break; // If an argument is in scope, prefer it
                    }
                }
                if (varindex == -1) {
                    printf("ERROR: variable %.*s not found!\n", subscript-start+1, latex+start);
                    exit(EXIT_FAILURE);
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
                    func_pos += parse_latex_rec(latex+subscript+7, arg1_end-subscript-13, function_list+func_pos, stack, variable_list, stack_size, &flags);
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
        if ((i+4 < end) && (strncmp(latex+i, "\\cdot", 5) == 0)) {
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

void parse_latex(char *latex, function *function_list, double *stack, variable *variable_list, int *stack_size) {
    // Find all of the variables in the expression
    int end = strlen(latex);
    int i=0;
    // Parse LaTeX
    uint8_t flags = 0;
    int func_pos = parse_latex_rec(latex, end, function_list, stack, variable_list, stack_size, &flags);

    for (i=0; i < func_pos; i++) {
        printf("%02d %p: first_arg: %p, next_arg: %p, oper: %p, value: %p\n", i, function_list+i, function_list[i].first_arg, function_list[i].next_arg, function_list[i].oper, function_list[i].value);
    }

    for (i=0; variable_list[i].name; i++)
        printf("%s variable pointer at %p, points to %p, type %08X\n", variable_list[i].name, &(variable_list[i].pointer), variable_list[i].pointer, variable_list[i].type);
}

expression* topological_sort_rec(expression *start) {
    expression *top = start;
    for (int i=0; i < start->num_dependencies; i++) {
        (start->dependencies)[i]->num_dependents--;
        if ((start->dependencies)[i]->num_dependents == 0) {
            // Recurse on a node
            (start->dependencies)[i]->next_expr = top;
            top = topological_sort_rec((start->dependencies)[i]);
        }
    }
    return top;
}

expression* topological_sort(expression *expression_list, int n_expr) {
    expression *top = NULL;
    for (int i=0; i < n_expr; i++) {
        if (expression_list[i].num_dependents == 0) {
            // Node may be used as a starting point
            expression_list[i].next_expr = top;
            top = topological_sort_rec(expression_list+i);
        }
    }
    return top;
}

int parse_file(function *function_list, double *stack, variable *variable_list, char *stringbuf, expression *expression_list, uint32_t *n_func, uint32_t *n_var, uint32_t *n_expr) {
    FILE *fp;
    char *line=NULL;
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
    printf("variable_list: %p\nexpression_list: %p\nstack: %p\nstringpos: %p\n", variable_list, expression_list, stack, stringpos);

    fp = fopen("/Users/matthias/Documents/plotter/test.txt", "r");
    if (fp == NULL) return 1;
    while ((read = getline(&line, &len, fp)) != -1) {
        printf("Found expression of length %zd\n", read);
        if (read <= 1) continue;
        printf("    %s", line);
        line[read-1] = 0;
        exprpos[0] = new_expression(NULL, NULL, NULL, 0, 0);
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
                else {
                    exprpos++;
                    continue;
                }
                // Check if the function has already been defined
                for (int j=0; variable_list[j].name; j++) {
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
                int end = extract_double(line, i+1, stack+stack_size);
                strncpy(stringpos, line, i);
                if (end+2 == read) {
                    printf("    Variable has value %f\n", stack[stack_size]);
                    varpos[0] = new_variable(stringpos, 1<<8, VARIABLE_IN_SCOPE, stack+stack_size);
                    stack_size++;
                } else {
                    varpos[0] = new_variable(stringpos, 0, VARIABLE_IN_SCOPE, NULL);
                }
                exprpos->var = varpos;
                varpos++;
                stringpos += i+1;
                i++;
            } else {
                exprpos++;
                continue;
            }
            // In either case (variable or function definition, we need to finish up the expression struct.
            // This involves extracting the variables used in the definition. However, we can't do this
            // until we parse all of the expressions. Thus, we store the index to the expression struct.
            exprpos->expr_begin = i;

        }
        exprpos++;
    }
    printf("loop done, %p, %p\n", variable_list->name, stringpos);
    for (i=0; variable_list[i].name; i++)
        printf("%s variable pointer at %p, points to %p, type %08X\n", variable_list[i].name, variable_list+i, variable_list[i].pointer, variable_list[i].type);
    fseek(fp, 0, SEEK_SET);

    *n_expr = exprpos - expression_list;
    printf("A total of %d expressions found\n", *n_expr);
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
    while ((read = getline(&line, &len, fp)) != -1) {
        if (read <= 1) continue;
        if (line[read-1] == '\n') line[read-1] = 0;
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
                        if (strncmp((((exprpos->dependencies)[k])->var)->name, line+j, subscript-j) == 0) {
                            printf("    Variable already in dependency table\n");
                            already_exists = 1;
                            break;
                        }
                    }
                    if (!already_exists) {
                        // If the dependency has not already been added, search for the variable in the
                        // expression table. Note that this will skip variables that are not defined by the
                        // user. For example, it will miss x and y as well as function arguments
                        for (int k=0; k < *n_expr; k++) {
                            if ((expression_list[k].var) && (strncmp((expression_list[k].var)->name, line+j, subscript-j) == 0)) {
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
                                if (strncmp(arg->name, line+j, subscript-j) == 0) {
                                    already_exists = 1;
                                    break;
                                }
                                arg++;
                            }
                        }
                        if (!already_exists && (line[j] != 'x') && (line[j] != 'y')) {
                            if (local_variable) {
                                printf("ERROR: Variable %.*s not defined\n", subscript-j, line+j);
                                exit(EXIT_FAILURE);
                            }
                            exprpos->flags |= EXPRESSION_PLOTTABLE;
                            exprpos->flags &= ~EXPRESSION_FIXED;
                            strncpy(stringpos, line+j, subscript-j);
                            *varpos = new_variable(stringpos, 1<<8, VARIABLE_IN_SCOPE | VARIABLE_XLIKE, NULL);
                            local_variable = varpos;
                            varpos++;
                            stringpos += subscript-j+1;
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
            func_pos += parse_latex_rec(line+i, read-i, function_list+func_pos, stack, variable_list, &stack_size, &flags);
            function_list[last_pos].value_type |= TYPE_ABSOLUTE_ADDR;
            function_list[last_pos].next_arg = (function*)((exprpos->var)+1);
            // Take the arguments out of the local scope after parsing
            first_arg = (exprpos->var)+1;
            while (first_arg->flags & VARIABLE_ARGUMENT) {
                first_arg->flags &= ~VARIABLE_IN_SCOPE;
                first_arg++;
            }
        } else if (!(exprpos->var) || (exprpos->var->pointer == 0)) {
            exprpos->func = function_list+func_pos;
            func_pos += parse_latex_rec(line+i, read-i, function_list+func_pos, stack, variable_list, &stack_size, &flags);
        }
        
        if (local_variable) local_variable->flags &= (~VARIABLE_IN_SCOPE);

        exprpos++;
    }

    // Connect all of the function calls to their definitions instead of their declarations
    // Connect implicit variables to their actual value (for example, replace u in w=10^{u} with x)
    for (i=0; i < func_pos; i++) {
        if (function_list[i].oper == func_user_defined) {
            function_list[i].value = ((variable*)(function_list[i].value))->pointer;
        }
        if (function_list[i].value_type & 0x40) {
            if (((variable*)(function_list[i].value))->flags & VARIABLE_XLIKE) function_list[i].value = (void*)variable_list;
            else if (((variable*)(function_list[i].value))->flags & VARIABLE_YLIKE) function_list[i].value = (void*)(variable_list+1);
        }
    }

    // Print the expression table
    for (i=0; expression_list+i < exprpos; i++)
        printf("expression pointer at %p, function %p, variable %p, flags %02x, begin %d, dep %p, num %d\n", expression_list+i, expression_list[i].func, expression_list[i].var, expression_list[i].flags, expression_list[i].expr_begin, expression_list[i].dependencies, expression_list[i].num_dependencies);
    
    // Print the function table
    for (i=0; i < func_pos; i++) {
        printf("%02d %p: first_arg: %p, next_arg: %p, oper: %p, value: %p, value_type: %08x\n", i, function_list+i, function_list[i].first_arg, function_list[i].next_arg, function_list[i].oper, function_list[i].value, function_list[i].value_type);
    }
    

    // Use topological sort to determine the order in which variable and function assignments must be evaluated.
    for (i=0; i < deptable_ofs; i++) 
        deptable[i]->num_dependents++;
    for (i=0; expression_list+i < exprpos; i++) {
        if (expression_list[i].num_dependencies > 0) {
            printf("%p\n", expression_list+i);
            for (int j=0; j < expression_list[i].num_dependencies; j++) printf("    %p\n", expression_list[i].dependencies[j]);
        }
    }
    for (i=0; expression_list+i < exprpos; i++) {
        printf("%d expressions depend on %p\n", expression_list[i].num_dependents, expression_list+i);
    }
    expression *top_expr = topological_sort(expression_list, *n_expr);
    printf("Top is %p\n", top_expr);
    expression *expr = top_expr;
    while (expr) {
        printf("'%p', ", expr);
        expr = expr->next_expr;
    }
    printf("\n");
    
    // Evaluate all variables in accordance with the topological ordering
    uint32_t type;
    expr = top_expr;
    while (expr) {
        if ((expr->func) && (expr->var) && !(expr->var->flags & VARIABLE_FUNCTION) && !(expr->flags & EXPRESSION_PLOTTABLE)) {
            printf("evaluating variable %s\n", expr->var->name);
            expr->var->pointer = stack+stack_size;
            type = (expr->func->oper(expr->func, stack + stack_size));
            printf("    result %f stored to %p, has type %08x\n", stack[stack_size], stack+stack_size, type);
            expr->var->type = type;
            stack_size += (type>>8);
        }
        expr = expr->next_expr;
    }

    // Evaluate all expressions that are not definitions and are not dependent on x or y
    for (expr=expression_list; expr < exprpos; expr++) {
        if ((expr->flags & EXPRESSION_FIXED) && !(expr->var)) {
            printf("evaluating expression %p\n", expr);
            type = (expr->func->oper(expr->func, stack + stack_size));
            if (type & TYPE_POINT) expr->flags |= EXPRESSION_PLOTTABLE;
            printf("    result %f stored to %p, has type %08x\n", stack[stack_size], stack+stack_size, type);
            expr->value = stack + stack_size;
            expr->value_type = type;
            stack_size += (type>>8);
        }
    }

    for (i=0; variable_list[i].name; i++)
        printf("%s variable pointer at %p, points to %p, type %08X, flags %02X\n", variable_list[i].name, variable_list+i, variable_list[i].pointer, variable_list[i].type, variable_list[i].flags);

    *n_func = func_pos;
    *n_var = varpos - variable_list;

    // Close the file
    fclose(fp);
    if (line) free(line);
    return 0;
}

