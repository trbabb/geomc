CC        = clang++
AR        = ar
PREFIX    = /usr/local
INCLUDES  = .
BUILD_STD = c++20
M_ARCH    = native
CFLAGS    = -g -march=$(M_ARCH) -std=$(BUILD_STD) -O3 -Wall \
	-fmessage-length=0 -Wno-unused-local-typedef \
	# -fsanitize=address -fno-omit-frame-pointer
LDFLAGS   = -g -lgeomc -lboost_unit_test_framework \
	# -fsanitize=address -fno-omit-frame-pointer
IFLAGS    = $(addprefix -I, $(INCLUDES))

LIB_SRC   = $(wildcard geomc/*.cpp) $(wildcard geomc/*/*.cpp)
TEST_SRC  = $(wildcard regression/*cpp)
ALL_SRC   = $(LIB_SRC) $(TEST_SRC) test/Profile.cpp

H_FILES   = $(shell find geomc -name "*.h")

TEST_BINS  = $(addprefix bin/, $(basename $(TEST_SRC)))
LIB_OBJS   = $(addprefix build/,    $(LIB_SRC:.cpp=.o))
EMLIB_OBJS = $(addprefix build/em/, $(LIB_SRC:.cpp=.o))

LIBNAME   = libgeomc.a
LIB       = lib/$(LIBNAME)
LITELIB   = lib/libgeomc_lite.a
EMLIB     = lib/libgeomc_em.a
INSTALL_INC_DIR = $(PREFIX)/include
INSTALL_LIB_DIR = $(PREFIX)/lib
WASM_INC_DIR = /usr/local/wasm/include
WASM_LIB_DIR = /usr/local/wasm/lib


all: lib litelib test

docs:
	mkdir -p doc/gen
	doxygen

clean:
	rm -rf ./build/*
	rm -f  ./lib/*.a
	rm -rf ./doc/gen/html

install: $(LIB) $(LITELIB)
	install -d $(foreach f, $(H_FILES), $(INSTALL_INC_DIR)/$(dir $(f)))
	install -d $(INSTALL_LIB_DIR)
	install -m 644 $(LIB) $(INSTALL_LIB_DIR)
	install -m 644 $(LITELIB) $(INSTALL_LIB_DIR)
	for f in $(H_FILES); do install -m 644 $$f $(INSTALL_INC_DIR)/$$f; done

install-wasm:
	mkdir -p $(WASM_INC_DIR)
	cp -rf ./geomc $(WASM_INC_DIR)
	cp -rf $(EMLIB) $(WASM_LIB_DIR)/$(LIBNAME)

uninstall:
	rm -rf $(INSTALL_INC_DIR)/geomc
	rm -f  $(INSTALL_LIB_DIR)/$(LIBNAME)
	rm -f  $(INSTALL_LIB_DIR)/$(notdir $(LITELIB))

test: $(TEST_BINS)
	$(foreach x, $(TEST_BINS), $(x);)

test-%: bin/regression/%
	$<

# basic library
lib: $(LIB)

# minimal library
litelib: $(LITELIB)

# emscripten-compatible library
emlib: $(EMLIB)


# auto-dependencies:
DEP = $(patsubst %.cpp, build/%.d, $(ALL_SRC))

-include $(DEP)


# build rules:

build/%.o: %.cpp build/%.d
	$(CC) -c $(CFLAGS) $(IFLAGS) -o $@ $<

build/em/%.o: %.cpp build/%.d
	@mkdir -p build/em/$(dir $<)
	emcc -c $(CFLAGS) $(IFLAGS) -o $@ $<

build/%.d: %.cpp
	@mkdir -p build/$(dir $<)
	$(CC) $(CFLAGS) $(IFLAGS) $< -MM -MT $(@:.d=.o) > $@

$(LIB): $(LIB_OBJS)
	@mkdir -p lib
	$(AR) rs $(LIB) $(LIB_OBJS)

$(LITELIB) : build/geomc/GeomException.o
	$(AR) rs $(LITELIB) build/geomc/GeomException.o

$(EMLIB) : $(EMLIB_OBJS)
	emar rs $(EMLIB) $(EMLIB_OBJS)

bin/regression/%: build/regression/%.o lib
	@mkdir -p $(dir $@)
	$(CC) $(LDFLAGS) $< -o $@

bin/%: build/%.o lib
	@mkdir -p $(dir $(patsubst build/%, bin/%, $<))
	$(CC) -g $< $(LIB) -o $@

