AUTOMAKE_OPTIONS = gnu

dist_doc_DATA = AUTHORS README

SUBDIRS = src

.PHONY: test	
test: check
	cd tests && $(MAKE) $(AM_MAKEFLAGS) test-local


