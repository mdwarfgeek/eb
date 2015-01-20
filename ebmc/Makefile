SRCS=main.c common.c constants.c fit.c func.c mc.c parms.c plots.c
OBJS=${SRCS:%.c=%.o}

include ../make.inc

# Compiler flags
CFLAGS+= -I../src -I../../lib $(PGPLOT_INC) $(LM_INC)

# Linker flags
LIBS=$(PGPLOT_LIBS) $(LM_LIBS) $(LAPACK_LIBS) -lm

all: ebmc

depend:
	$(CC) $(CFLAGS) -E -MM $(SRCS) > .depend

../src/libeb.a:
	(cd ../src && $(MAKE) libeb.a)

../../lib/liblfa.a:
	(cd ../../lib && $(MAKE))

ebmc: $(OBJS) ../src/libeb.a ../../lib/liblfa.a
	$(FC) -o $@ $(OBJS) ../src/libeb.a ../../lib/liblfa.a $(LIBS)

clean:
	rm -f $(OBJS) ebmc
	rm -f .depend
