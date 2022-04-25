PYTHON?=python
LIB?=eb*.so

.PHONY: all test clean

all: $(LIB)

../src/libeb.a:
	(cd ../src && $(MAKE))

$(LIB): wrap.c ../src/libeb.a
	$(PYTHON) setup.py build_ext --inplace

test: $(LIB)
	$(PYTHON) test.py

clean:
	rm -f $(LIB)
	rm -rf build
