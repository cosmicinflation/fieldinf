AUTOMAKE_OPTIONS = gnu

ACLOCAL_AMFLAGS= -I m4

dist_doc_DATA = AUTHORS README

SUBDIRS = src

.PHONY: test	
test: check
	cd tests && $(MAKE) $(AM_MAKEFLAGS) test-local


