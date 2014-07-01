all:
	(cd src && $(MAKE))
	(cd python && $(MAKE))

depend:
	(cd src && $(MAKE) depend)

clean:
	(cd src && $(MAKE) clean)
	(cd python && $(MAKE) clean)
