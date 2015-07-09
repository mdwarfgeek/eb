# Makefile for EB model library

C_SRCS=eb_model_flt.c eb_model_dbl.c eb_strings.c eb_util.c \
	dellke_gen.c dellpi_gen.c \
	fellke_gen.c fellpi_gen.c
S_SRCS=
OBJS=${C_SRCS:%.c=%.o} ${S_SRCS:%.S=%.o}

include ../make.inc

all: machconst.h libeb.a

depend: machconst.h
	$(CC) $(CFLAGS) -E -MM $(C_SRCS) $(S_SRCS) > .depend

machconst.h: makemach
	./makemach > machconst.h

makemach: makemach.c
	$(CC) $(CFLAGS) -o makemach makemach.c -lm

libeb.a: $(OBJS)
	rm -f $@
	ar r $@ $(OBJS)
	ranlib $@

clean:
	rm -f $(OBJS) libeb.a
	rm -f machconst.h makemach
	rm -f .depend
