AM_CFLAGS = -I/$(top_srcdir)/src  $(GMP_CFLAGS) $(ATLAS_CFLAGS)

TESTLIB = $(top_builddir)/src/libiml.la $(ATLAS_LIBS) $(GMP_LIBS)

TESTS = test-smallentry$(EXEEXT) test-largeentry$(EXEEXT)

EXTRA_PROGRAMS = test-smallentry test-largeentry


test_smallentry_SOURCES = test-smallentry.c 
test_smallentry_CFLAGS = $(AM_CFLAGS)
test_smallentry_LDADD = $(TESTLIB)


test_largeentry_SOURCES = test-largeentry.c 
test_largeentry_CFLAGS = $(AM_CFLAGS)
test_largeentry_LDADD = $(TESTLIB)


