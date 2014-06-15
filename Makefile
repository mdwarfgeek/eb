all:
	(cd src && $(MAKE))
	(cd python && $(MAKE))

clean:
	(cd src && $(MAKE) clean)
	(cd python && $(MAKE) clean)
