################################################################################
# Makefile for despyfits
################################################################################

################################################################################
# Install and include locations.
################################################################################

SHELL=/bin/sh
INCLUDE := ./include
LIB := ./lib
SRC := ./src

################################################################################
# Headers changes in which should force rebuilds
################################################################################
HEADERS := despyfits/mask_bits.h despyfits/desimage.h

################################################################################
# Explore our environment
################################################################################

# Define Flavor Shared lib extension
FLAVOR = $(shell uname -s) # this give either Darwin/Linux
$(info Discovered FLAVOR: $(FLAVOR))
ifneq (,$(findstring Darwin, ${FLAVOR}))
   SHLIB_SUFFIX = dylib
endif 
ifneq (,$(findstring Linux, ${FLAVOR}))
   SHLIB_SUFFIX = so
endif 

# Paths to search for dependecies:
vpath %.h $(INCLUDE)
vpath %.o $(LIB)

################################################################################
# C flags
################################################################################
CC :=		gcc
CFLAGS := 	-O3 -g -Wall -shared -fPIC

override CFLAGS += $(addprefix -I ,$(INCLUDE))
override CFLAGS += $(addprefix -L ,$(LIBS))

################################################################################
# The libraries we are to make
################################################################################
SHLIBS = $(LIB)/libbpm.$(SHLIB_SUFFIX) \
	 $(LIB)/libfixcol.$(SHLIB_SUFFIX) \
	 $(LIB)/libmasksatr.$(SHLIB_SUFFIX) \
	 $(LIB)/libfoo.$(SHLIB_SUFFIX) \
	 $(LIB)/libbar.$(SHLIB_SUFFIX)

################################################################################
# The rule to actually make the libraries
################################################################################
$(LIB)/%.$(SHLIB_SUFFIX) : $(SRC)/%.c $(HEADERS)
	$(CC) $(CFLAGS) $< -o $@

################################################################################
# Generic rules
################################################################################
.PHONY: all
all: $(SHLIBS)

tidy:
	rm -rf *.o *~ \#*\#

clean: tidy
	for f in $(SHLIBS) ; do rm -f $$f ; done

install:
	false #nothing to install here.
