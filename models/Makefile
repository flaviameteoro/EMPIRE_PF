SUBDIRS := $(wildcard */.)

all: models

models:
	-for d in $(SUBDIRS); do (cd $$d; $(MAKE) -e ); done
