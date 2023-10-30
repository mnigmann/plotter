#include <gtk/gtk.h>
#include <stdio.h>
#include <stdint.h>
#include "parse.h"
#include "functions.h"

function *function_list;
function *orig_addr;
uint32_t n_func;

enum
{
  COL_POINTER = 0,
  COL_OPER,
  COL_VALUE,
  COL_VALUE_TYPE,
  COL_VARIABLE,
  NUM_COLS
} ;

void insert_rec(GtkTreeStore *store, function *func, GtkTreeIter *parent, char *stringbuf) {
    GtkTreeIter iter;
    printf("inserting %p, %p, %p to %p, %p, %p\n", func, func->oper, func->value, store, &iter, parent);
    uint32_t v1 = sprintf(stringbuf, "%p", func)+1;
    oper_data *oper = oper_lookup(func->oper);
    uint32_t v2;
    if (oper) v2 = v1+sprintf(stringbuf+v1, "%s", oper->name)+1;
    else v2 = v1+sprintf(stringbuf+v1, "%p", func->oper)+1;
    uint32_t v3 = v2+sprintf(stringbuf+v2, "%p", func->value)+1;
    uint32_t v4 = v3+sprintf(stringbuf+v3, "%08x", func->value_type)+1;
    stringbuf[v4] = 0;
    if (func->value_type & 0x40) sprintf(stringbuf+v4, "%s", ((variable*)(func->value))->name);
    if (func->oper == func_user_defined) sprintf(stringbuf+v4, "%s", (((variable*)(((function*)(func->value))->next_arg))-1)->name);
    gtk_tree_store_append(store, &iter, parent);
    gtk_tree_store_set(store, &iter, COL_POINTER, stringbuf, COL_OPER, stringbuf+v1, COL_VALUE, stringbuf+v2, COL_VALUE_TYPE, stringbuf+v3, COL_VARIABLE, stringbuf+v4, -1);
    function *arg = func->first_arg;
    while (arg) {
        insert_rec(store, arg, &iter, stringbuf);
        arg = arg->next_arg;
    }
}

static GtkTreeModel *create_and_fill_model(file_data *fd) {
    printf("activating 3\n");
    GtkTreeStore *store = gtk_tree_store_new(NUM_COLS, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING, G_TYPE_STRING);

    /* Append a row and fill in some data */
    GtkTreeIter iter;
    char stringbuf[200];
    
    expression *expr;
    for (int i=0; i < fd->n_expr; i++) {
        expr = (fd->expression_list)+i;
        uint32_t v1 = sprintf(stringbuf, "Expression %p", expr)+1;
        uint32_t v2 = v1+1; stringbuf[v1] = 0;
        uint32_t v3 = v2+sprintf(stringbuf+v2, "%p", expr->value)+1;
        uint32_t v4 = v3+sprintf(stringbuf+v3, "%08x", expr->value_type)+1;
        stringbuf[v4] = 0;
        if (expr->var) sprintf(stringbuf+v4, "%s", expr->var->name);
        gtk_tree_store_append(store, &iter, NULL);
        gtk_tree_store_set(store, &iter, COL_POINTER, stringbuf, COL_OPER, stringbuf+v1, COL_VALUE, stringbuf+v2, COL_VALUE_TYPE, stringbuf+v3, COL_VARIABLE, stringbuf+v4, -1);
        insert_rec(store, expr->func, &iter, stringbuf);
    }

    return GTK_TREE_MODEL(store);
}

void treeview_activate(file_data *fd) {
    printf("activating 1\n");
    GtkWidget *window;

    window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    GtkWidget *scroll = gtk_scrolled_window_new(NULL, NULL);
    GtkWidget *treeview = gtk_tree_view_new();
    gtk_container_add(GTK_CONTAINER(scroll), treeview);
    gtk_container_add(GTK_CONTAINER(window), scroll);
    GtkCellRenderer *renderer;
    printf("activating 2\n");
    renderer = gtk_cell_renderer_text_new ();
    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW (treeview), -1, "Pointer", renderer,
                                                "text", COL_POINTER,
                                                NULL);

    renderer = gtk_cell_renderer_text_new ();
    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW (treeview), -1, "Operation", renderer,
                                                "text", COL_OPER,
                                                NULL);

    renderer = gtk_cell_renderer_text_new ();
    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW (treeview), -1, "Value", renderer,
                                                "text", COL_VALUE,
                                                NULL);

    renderer = gtk_cell_renderer_text_new ();
    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW (treeview), -1, "Value type", renderer,
                                                "text", COL_VALUE_TYPE,
                                                NULL);

    renderer = gtk_cell_renderer_text_new ();
    gtk_tree_view_insert_column_with_attributes(GTK_TREE_VIEW (treeview), -1, "Variable", renderer,
                                                "text", COL_VARIABLE,
                                                NULL);

    GtkTreeModel *model = create_and_fill_model(fd);

    gtk_tree_view_set_model(GTK_TREE_VIEW(treeview), model);

    /* The tree view has acquired its own reference to the
     *  model, so we can drop ours. That way the model will
     *  be freed automatically when the tree view is destroyed
     */
    g_object_unref(model);

    gtk_widget_show_all(window);
}

/*
int main(int argc, char **argv) {
    FILE *fp = fopen("/tmp/function_list", "r");
    fseek(fp, 0, SEEK_END);
    long length = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    fread(&orig_addr, sizeof(int64_t), 1, fp);
    length -= sizeof(int64_t);
    n_func = length / sizeof(function);
    printf("loaded %d functions\n", n_func);
    function_list = malloc(length);
    fread(function_list, 1, length, fp);
    fclose(fp);
    // Print the function table
    int i;
    for (i=0; i < 10; i++) {
        printf("%02d %p: first_arg: %p, next_arg: %p, oper: %p, value: %p, value_type: %08x\n", i, function_list+i, function_list[i].first_arg, 
            function_list[i].next_arg, function_list[i].oper, function_list[i].value, function_list[i].value_type);
    }

    GtkApplication *app;
    int status;


    app = gtk_application_new("org.gtk.example", G_APPLICATION_DEFAULT_FLAGS);
    g_signal_connect(app, "activate", G_CALLBACK (activate), NULL);
    status = g_application_run(G_APPLICATION (app), argc, argv);
    g_object_unref (app);

    return status;
}*/
