CC = clang++
AR = ar
PREFIX   = /opt/local
INCLUDES = .
CFLAGS   = -std=c++11 -O3 -Wall -c -fmessage-length=0 
IFLAGS   = $(addprefix -I, $(INCLUDES))
MODULES  = function linalg random shape
SOURCES  = $(wildcard geomc/*.cpp) \
           $(foreach m, $(MODULES), $(wildcard geomc/$(m)/*.cpp))

OBJECTS  = $(addprefix build/, $(notdir $(SOURCES:.cpp=.o)))
LIBNAME  = libgeomc.a
LITELIB  = libgeomc_lite.a
LIB      = lib/$(LIBNAME)
INCDIR   = $(PREFIX)/include
LIBDIR   = $(PREFIX)/lib

all : lib liblite

docs :
	mkdir -p doc/gen
	doxygen

lib : $(OBJECTS)
	mkdir -p lib
	$(AR) rs $(LIB) $(OBJECTS)
	@echo
	@echo Done building library.

liblite : build/GeomException.o build/Hash.o
	$(AR) rs lib/$(LITELIB) build/GeomException.o build/Hash.o

profile : lib build/Profile.o
	mkdir -p bin
	$(CC) -g build/Profile.o $(LIB) -o bin/profile

build/Profile.o : test/Profile.cpp
	$(CC) -g $(CFLAGS) $(IFLAGS) -o build/Profile.o test/Profile.cpp 

build/%.o : geomc/random/%.cpp
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $<

build/%.o : geomc/%.cpp
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $<

install : all
	mkdir -p $(INCDIR)
	cp -rf ./geomc $(INCDIR)
	cp -rf $(LIB) $(LIBDIR)
	cp -rf lib/$(LITELIB) $(LIBDIR)

clean :
	rm -f  ./build/*.o
	rm -f  ./lib/*.a
	rm -rf ./doc/gen/html

uninstall :
	rm -rf $(INCDIR)/geomc
	rm -f $(LIBDIR)/$(LIBNAME)
	rm -f $(LIBDIR)/$(LITELIB)

