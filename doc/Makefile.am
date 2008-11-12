## Process this file with automake to generate Makefile.in
##

EXTRA_DIST = docgen.in document.pl Changes.html Doxyfile.in style.css GettingStarted.html

install-data-local: 
	$(INSTALL) -d html
	sed -e 's%ZzZzZ%$(VERSION)%' docgen.in > docgen.pl
	perl docgen.pl  html
	cp style.css     html/.
	cp Changes.html html/.
	cp GettingStarted.html html/.
	sed -e 's%ZzZzZ%html/api%' Doxyfile.in > Doxyfile
	doxygen Doxyfile || :
	chmod 0755 html
	@for file in `find html`; do \
	if test -d $$file; then \
	  chmod 0755 $$file; \
	else \
	  chmod 0644 $$file; \
	fi; \
	done
	-rm -r     $(prefix)/html/$(PACKAGE)
	-mkdir -p  $(prefix)/html
	cp -r html $(prefix)/html/$(PACKAGE)

uninstall-local:
	rm -r $(prefix)/html/$(PACKAGE)

clean-local:
	-rm -r html