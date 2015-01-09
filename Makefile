SHELL=/bin/sh

DIRS=src test

#
# The pattern of this makefile is to call a subordinate makefile for
# each of the tasks, if a Makefile is is present in the subordinate
# directory, else to to the task. Install is an exception.
#

all: build localinstall

build: 
	for d in $(DIRS) ; do if [ -e $$d/Makefile ] ; then (cd $$d;$(MAKE) $(MAKEFLAGS) CFLAGS="$(CFLAGS)" CC="$(CC)" all ) fi ; done

clean: tidy
	for d in $(DIRS) ; do if [ -e $$d/Makefile ] ; then (cd $$d;$(MAKE) $(MAKEFLAGS) clean  ) fi ; done

tidy:
	for d in $(DIRS) ; do if [   -e $$d/Makefile ] ; then (cd $$d;$(MAKE) $(MAKEFLAGS) tidy  )  fi ; done
	for d in $(DIRS) ; do if [ ! -e $$d/Makefile ] ; then (cd $$d ; rm -f \#*\#  *~ )           fi ; done
	rm -rf \#*\#  *~

localinstall:
	cp src/*.so lib
