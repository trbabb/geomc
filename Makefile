CC = g++
AR = ar
INCLUDES = /usr/local/boost .
CFLAGS   = -O3 -Wall -c -fmessage-length=0 -Wno-unused -Wno-unused-local-typedefs
IFLAGS   = $(addprefix -I, $(INCLUDES))
MODULES  = function linalg random shape
SOURCES  = $(wildcard geomc/*.cpp) \
           $(foreach m, $(MODULES), $(wildcard geomc/$(m)/*.cpp))

OBJECTS  = $(addprefix build/, $(notdir $(SOURCES:.cpp=.o)))
LIBNAME  = libgeomc.a
LIB      = lib/$(LIBNAME)
INCDIR   = /opt/local/include
LIBDIR   = /opt/local/lib

all : lib

docs :
	mkdir -p doc/gen
	doxygen

lib : $(OBJECTS)
	$(AR) rs $(LIB) $(OBJECTS)
	@echo
	@echo Done building library.

profile : lib build/Profile.o
	$(CC) $(LIB) build/Profile.o -o bin/profile

build/Profile.o : test/Profile.cpp
	$(CC) $(CFLAGS) $(IFLAGS) -o build/Profile.o test/Profile.cpp 

build/%.o : geomc/random/%.cpp
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $<

build/%.o : geomc/%.cpp
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $<

install : all
	mkdir -p $(INCDIR)
	cp -rf ./geomc $(INCDIR)
	cp -rf $(LIB) $(LIBDIR)

clean :
	rm -f  ./build/*.o
	rm -f  ./lib/*.a
	rm -rf ./doc/gen/html

uninstall :
	rm -rf $(INCDIR)/geomc
	rm -f $(LIBDIR)/$(LIBNAME)
