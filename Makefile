BUILD_DIR := ./build
SRC_DIR := .

SRCS := $(wildcard *.c)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

LDFLAGS := `pkg-config --libs gtk+-3.0` -lm -L/usr/local/opt/lapack/lib -llapack -L/usr/local/lib -lgsl -lgslcblas
GCCFLAGS := `pkg-config --cflags gtk+-3.0` -I/usr/local/opt/lapack/include -I/usr/local/include

$(BUILD_DIR)/plotter: $(OBJS)
	gcc $(OBJS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.c.o: %.c
	gcc $(GCCFLAGS) -c $< -o $@

