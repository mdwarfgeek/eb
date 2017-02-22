IDL?=$(IDL_DIR)/bin/idl
LIB?=eb.so

.PHONY: all clean

all: $(LIB)

../src/libeb.a:
	(cd ../src && $(MAKE))

$(LIB): wrap.c ../src/libeb.a
	echo "make_dll,['wrap'],'eb','IDL_Load',input_directory='./',extra_cflags='-I../src -Wall',extra_lflags='../src/libeb.a -lm',compile_directory='./',/show_all_output,output_directory='./'" | $(IDL)

clean:
	rm -f $(LIB)
