# Makefile for EB model library

C_SRCS=eb_model_flt.c eb_model_dbl.c eb_strings.c eb_util.c \
	dellke_gen.c dellppi_gen.c fellke_gen.c fellppi_gen.c
S_SRCS=
OBJS=${C_SRCS:%.c=%.o} ${S_SRCS:%.S=%.o}

include ../make.inc

all: libeb.a

depend:
	$(CC) -E -MM $(C_SRCS) $(S_SRCS) > .depend

libeb.a: $(OBJS)
	rm -f $@
	ar r $@ $(OBJS)
	ranlib $@

clean:
	rm -f $(OBJS) libeb.a
	rm -f .depend

sinclude .depend
