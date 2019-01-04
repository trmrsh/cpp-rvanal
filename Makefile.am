## Process this file with automake to generate Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src doc include

## file of aliases
ALIASES = Rvanal

.PHONY : $(ALIASES)

DATE    = $(shell date)

progdir = @bindir@/@PACKAGE@

prog_DATA = $(ALIASES)

# Note horribly complex looking line to generate a file that can be sourced
# unde both bash and csh. It's designed to get through variable expansion
$(ALIASES):
	@echo '#' > $(ALIASES)
	@echo '# This file defines aliases for the $(PACKAGE) v$(VERSION) software.' >> $(ALIASES)
	@echo '# Source it to define them.' >> $(ALIASES)
	@echo '#' >> $(ALIASES)
	@echo '# this to allow both bash and csh to work' >> $(ALIASES)
	@echo 'test "$$?BASH_VERSION" = "0" || eval '\''alias() { command alias "$$1=$$2"; }'\' >> $(ALIASES)
	@echo '#' >> $(ALIASES)
	@echo 'alias aliases   $(progdir)/aliases'   >> $(ALIASES)
	@echo 'alias bayes     $(progdir)/bayes'     >> $(ALIASES)
	@echo 'alias dprob     $(progdir)/dprob'     >> $(ALIASES)
	@echo 'alias mergers   $(progdir)/mergers'   >> $(ALIASES)
	@echo 'alias nthalias  $(progdir)/nthalias'  >> $(ALIASES)
	@echo 'alias pdist     $(progdir)/pdist'     >> $(ALIASES)
	@echo 'alias pfold     $(progdir)/pfold'     >> $(ALIASES)
	@echo 'alias pgram     $(progdir)/pgram'     >> $(ALIASES)
	@echo 'alias sindet    $(progdir)/sindet'    >> $(ALIASES)
	@echo 'alias sinfalse  $(progdir)/sinfalse'  >> $(ALIASES)
	@echo 'alias sintest   $(progdir)/sintest'   >> $(ALIASES)
	@echo 'alias vtest     $(progdir)/vtest'     >> $(ALIASES)
	@echo ' ' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Welcome to $(PACKAGE) version $(VERSION), built $(DATE)."' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Commands available are: "' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "aliases, bayes, dprob, mergers, nthalias, pdist, pfold, pgram, sindet, sinfalse, sintest and vtest"' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Point a browser at ${prefix}/html/$(PACKAGE)/index.html for help."' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
