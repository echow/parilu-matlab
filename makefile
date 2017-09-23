# make sure this includes the architectures you care about
SYS = $(shell uname -s)
ARCH = $(shell uname -m)
ifeq ($(SYS), Linux)
  MEX = mex
  ifeq ($(ARCH), x86_64)
    MEXSUFFIX = mexa64
  endif
  ifeq ($(ARCH), i686)
    MEXSUFFIX = mexglx
  endif
else # CYGWIN
  MEX = mex.bat
  MEXSUFFIX = mexw64
endif

# set sources explicitly if not all your .c files are mex programs
SOURCES = $(wildcard *.c)
TARGETS = $(SOURCES:.c=.$(MEXSUFFIX))

all: $(TARGETS)

%.$(MEXSUFFIX): %.c
	$(MEX) -O -largeArrayDims CFLAGS='$$CFLAGS -fopenmp -std=c99' LDFLAGS='$$LDFLAGS -fopenmp' $<

clean:
	rm -f *.$(MEXSUFFIX)

.PHONY: clean
