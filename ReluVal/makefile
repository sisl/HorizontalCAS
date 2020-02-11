CC=gcc
CFLAGS=-I /opt/OpenBLAS/include/ 
LDFLAGS= -L/opt/OpenBLAS/lib -lopenblas -lpthread -lm
DEPS=matrix.h nnet.h split.h

all: network_test_HCAS
all: CFLAGS += -O3
all: CFLAGS += -std=c99
all: LDFLAGS += -O3

debug: network_test_HCAS
debug: CFLAGS += -DDEBUG -g
debug: CFLAGS += -std=c99
debug: LDFLAGS += -DDEBUG -g

network_test_HCAS: matrix.o nnet.o network_test_HCAS.o split.o
	$(CC) -static -o $@ $^ $(LDFLAGS)

c.o: 
	$(CC) $(CFLAGS) $<  -o $@

clean:
	rm -f *.o network_test_HCAS

