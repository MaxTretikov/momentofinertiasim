SRCS=voxelize.c voxellib.c
DEPS=voxellib.h
OBJS=$(SRCS:.c=.o)
TARGET=voxelize
CFLAGS=-DISPYTHON=0 -g -Wall -O0 -fPIC
LIBS=-lm

all: $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(LIBS)
