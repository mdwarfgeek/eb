.PHONY: all python idl depend clean

all:
	(cd src && $(MAKE))
	(cd python && $(MAKE))

python:
	(cd python && $(MAKE))

idl:
	(cd idl && $(MAKE))

depend:
	(cd src && $(MAKE) depend)

clean:
	(cd src && $(MAKE) clean)
	(cd python && $(MAKE) clean)
	(cd idl && $(MAKE) clean)
