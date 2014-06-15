PYTHON?=python
LIB?=eb.so

all: $(LIB)

../src/libeb.a:
	(cd ../src && $(MAKE) libeb.a)

$(LIB): wrap.c ../src/libeb.a
	$(PYTHON) setup.py build_ext --inplace

test: $(LIB)
	$(PYTHON) test.py

clean:
	rm -f $(LIB)
	rm -rf build
