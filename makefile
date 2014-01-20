DIRS    := src
SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.c)) src/gitversion.c
OBJS    := $(patsubst %.c, %.o, $(SOURCES))
OBJS    := $(foreach o,$(OBJS),./obj/$(o))
DEPFILES:= $(patsubst %.o, %.P, $(OBJS))
 
CFLAGS   = -Wall -MMD -c -O3 -std=c99
INCLUDES = -I/usr/include/libxml2 -I./inc
LIBS     = -lm -lxml2 -lpthread
CC       = gcc

all: pepquant2

#link the executable
pepquant2: $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o pepquant2
	rm -f src/gitversion.c
 
#generate dependency information and compile
obj/%.o : %.c 
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -MF obj/$*.P $<
	@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < obj/$*.P >> obj/$*.P;

src/gitversion.c: .git/HEAD .git/index
	echo "const char *gitversion = \"$(shell git rev-parse HEAD)\";" > $@
	echo "const char *commit = \"$(shell git rev-list HEAD --count)\";" >> $@
 
#remove all generated files
clean:
	rm -f main
	rm -rf obj
	rm -f src/gitversion.c

 
#include the dependency information
-include $(DEPFILES)

#makefile outline based on
#http://andreaslindh.wordpress.com/2011/
#10/05/auto-dependency-makefile-example/
