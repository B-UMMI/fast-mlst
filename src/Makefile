######### IMPORTANT DO NOT COMMENT THESE FLAGS

# The compiler to use
CC        = gcc
# Ansi flags, make the compiler more informative
STANDARD  = -Wall -Wextra -std=c99 -I../divsufsort
# Puts cacheline size in a macro
CLS       = -DCLS=$(getconf LEVEL1_DCACHE_LINESIZE)
# Use this if you want to use SSE2 instructions
SSE2      = -msse2
# Debug flags for gdb/GUD
DEBUG     = -O0 -gdwarf-2 -g3 -ggdb
# Profile flags for gprof
PROFILE   = -pg -fprofile-arcs -ftest-coverage
# Optimization flags for production
OPTIMIZE  = -O3 -DNDEBUG
# Static for static StandAlones
STATIC = -static
# PThreads are usefull
THREADS = -pthread

#########  COMMENT THESE

# Uncomment the flags that you need in a particular compile
CFLAGS   = $(STANDARD)
#CFLAGS  += $(DEBUG)
#CFLAGS  += -D_BSD_SOURCE
#CFLAGS  += $(PROFILE)
CFLAGS  += $(OPTIMIZE)
#CFLAGS  += $(CLS) $(SSE2)
#CFLAGS  += $(STATIC)
#CFLAGS  += $(THREADS)

LDFLAGS = -lm

# Files
EXECS = main

main_CS = main.c qsufsort.c sautils.c
main_HS = main.h sautils.h
main_OS = main.o qsufsort.o sautils.o

# Phony targets 
.PHONY: all clean depend

# Default Compile
all: $(EXECS)

## Linking rules
main: $(main_OS)
	@echo Linking: $@
	$(CC) $(CFLAGS) $(main_OS) -o $@ $(LDFLAGS)

## Build Object
%.o: %.c
	@echo Build Object from: $<
	$(CC) $(CFLAGS) -c -o $@ $<

depend.mak: $(test_CS) $(test_HS) $(main_CS) $(test_HS)
	@echo Making dependencies ...
	$(CC) $(CFLAGS) -MM $(test_CS) $(main_CS) > depend.mak

-include depend.mak

## Clean up
clean:
	@echo Cleaning Up
	rm -f $(EXECS) *.o gmon.out depend.mak *~ callgrind.out.* cachegrind.out.*

