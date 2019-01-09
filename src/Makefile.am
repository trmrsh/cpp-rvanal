## Process this file with automake to generate Makefile.in
##
## This is the file that must be edited if you are changing anything in the source directory

## Programs

progdir = @bindir@/@PACKAGE@

prog_PROGRAMS     = aliases bayes dprob nthalias pfold pdist pgram sindet sinfalse sintest vtest mergers

aliases_SOURCES  = aliases.cc
bayes_SOURCES    = bayes.cc
dprob_SOURCES    = dprob.cc
nthalias_SOURCES = nthalias.cc
mergers_SOURCES  = mergers.cc
pdist_SOURCES    = pdist.cc
pfold_SOURCES    = pfold.cc
pgram_SOURCES    = pgram.cc
sindet_SOURCES   = sindet.cc
sinfalse_SOURCES = sinfalse.cc
sintest_SOURCES  = sintest.cc
vtest_SOURCES    = vtest.cc

AM_CPPFLAGS = -I../include -I../. $(PCRECPP_CPPFLAGS)

LDADD = librvanal.la $(PCRECPP_LIBS)

## Library

lib_LTLIBRARIES = librvanal.la 

librvanal_la_SOURCES = bayes_prob.cc bayes_trap.cc bestsin.cc dpb.cc sinfit_chisq.cc m2min.cc fm.cc \
sinfit_errors.cc

DATE = $(shell date)

install-data-hook:
	echo "This is $(PACKAGE)-$(VERSION), built on $(DATE)" > $(progdir)/VERSION




