## Process this file with automake to generate Makefile.in

SUBDIRS = src doc include

EXTRA_DIST = pgplot_libs

export: dist
	scp $(PACKAGE)-$(VERSION).tar.gz $(WEB_SERVE):$(WEB_PATH)/software/.
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software ; \rm -f $(PACKAGE).tar.gz ; ln -s $(PACKAGE)-$(VERSION).tar.gz $(PACKAGE).tar.gz; \rm -fr $(PACKAGE); mkdir $(PACKAGE)"
	scp -r doc/html $(WEB_SERVE):$(WEB_PATH)/software/$(PACKAGE)/.
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software; sed -e 's%<\!-- rvanal version --><td>.*</td>%<\!-- rvanal version --><td>$(VERSION)</td>%' index.html > junk.html; mv junk.html index.html"
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software; sed -e 's%<\!-- rvanal date --><td>.*</td>%<\!-- rvanal date --><td>$(DATE)</td>%' index.html > junk.html; mv junk.html index.html"

## file of aliases

ALIASES = Rvanal

.PHONY : $(ALIASES)

DATE    = $(shell date)

progdir = @bindir@/@PACKAGE@

prog_DATA = $(ALIASES)

$(ALIASES):
	@echo '#' > $(ALIASES)
	@echo '# This file defines aliases for the $(PACKAGE) v$(VERSION) software.' >> $(ALIASES)
	@echo '# Source it to define them.' >> $(ALIASES)
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
