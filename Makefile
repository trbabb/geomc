CC        = g++
AR        = ar
PREFIX    = /opt/local
INCLUDES  = .
CFLAGS    = -g -std=c++11 -O3 -Wall -fmessage-length=0 -Wno-unused-local-typedef
IFLAGS    = $(addprefix -I, $(INCLUDES))

LIB_SRC   = $(wildcard geomc/*.cpp) $(wildcard geomc/*/*.cpp)
TEST_SRC  = $(wildcard regression/*cpp)
ALL_SRC   = $(LIB_SRC) $(TEST_SRC) test/Profile.cpp

TEST_BINS = $(addprefix bin/, $(basename $(TEST_SRC)))
LIB_OBJS  = $(addprefix build/, $(LIB_SRC:.cpp=.o))

LIBNAME   = libgeomc.a
LIB       = lib/$(LIBNAME)
LITELIB   = geomc_lite.a
INSTALL_INC_DIR = $(PREFIX)/include
INSTALL_LIB_DIR = $(PREFIX)/lib


all: lib liblite test

docs:
	mkdir -p doc/gen
	doxygen

clean:
	rm -rf ./build/*
	rm -f  ./lib/*.a
	rm -rf ./doc/gen/html

install: all
	mkdir -p $(INSTALL_INC_DIR)
	cp -rf ./geomc $(INSTALL_INC_DIR)
	cp -rf $(LIB) $(INSTALL_INC_DIR)
	cp -rf lib/$(LITELIB) $(INSTALL_LIB_DIR)

uninstall:
	rm -rf $(INSTALL_INC_DIR)/geomc
	rm -f  $(INSTALL_LIB_DIR)/$(LIBNAME)
	rm -f  $(INSTALL_LIB_DIR)/$(LITELIB)

test: $(TEST_BINS)
	$(foreach x, $(TEST_BINS), $(x);)

test-%: bin/regression/%
	$<


# auto-dependencies:
DEP = $(patsubst %.cpp, build/%.d, $(ALL_SRC))

-include $(DEP)


# build rules:

build/%.o: %.cpp build/%.d
	$(CC) -c $(CFLAGS) $(IFLAGS) -o $@ $<

build/%.d: %.cpp
	@mkdir -p build/$(dir $<)
	$(CC) $(CFLAGS) $(IFLAGS) $< -MM -MT $(@:.d=.o) > $@

lib: $(LIB_OBJS)
	@mkdir -p lib
	$(AR) rs $(LIB) $(LIB_OBJS)

liblite : build/geomc/GeomException.o build/geomc/Hash.o
	$(AR) rs lib/$(LITELIB) build/geomc/GeomException.o build/geomc/Hash.o

bin/regression/%: build/regression/%.o lib
	@mkdir -p $(dir $@)
	$(CC) -g -lgeomc -lboost_unit_test_framework $< -o $@

bin/%: build/%.o lib
	@mkdir -p $(dir $(patsubst build/%, bin/%, $<))
	$(CC) -g $< $(LIB) -o $@

