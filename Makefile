CC = g++
AR = ar
INCLUDES = /usr/local/boost .
CFLAGS   = -O3 -Wall -c -fmessage-length=0 -Wno-unused -Wno-unused-local-typedefs
IFLAGS   = $(addprefix -I, $(INCLUDES))
SOURCES  = geomc/GeomException.cpp \
	   geomc/random/LCRand.cpp \
	   geomc/random/MTRand.cpp \
	   geomc/random/Random.cpp \
	   geomc/random/RandomTools.cpp

OBJECTS  = $(addprefix build/, $(notdir $(SOURCES:.cpp=.o)))
LIB      = lib/libgeomc.a
INCDIR   = /opt/local/include
LIBDIR   = /opt/local/lib

all: lib
	echo done.

lib: $(OBJECTS)
	$(AR) rs $(LIB) $(OBJECTS)

profile: lib build/Profile.o
	$(CC) $(LIB) build/Profile.o -o bin/profile

build/Profile.o: test/Profile.cpp
	$(CC) $(CFLAGS) $(IFLAGS) -o build/Profile.o test/Profile.cpp 

build/%.o : geomc/random/%.cpp
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $<

build/%.o : geomc/%.cpp
	$(CC) $(CFLAGS) $(IFLAGS) -o $@ $<

install: all
	mkdir -p $(INCDIR)
	cp -rf ./geomc $(INCDIR)
	cp -rf $(LIB) $(LIBDIR)

clean:
	rm -f ./build/*.o
	rm -f ./lib/*.a

uninstall:
	rm -rf $(INCDIR)/geomc
	rm -f $(LIBDIR)/$(LIB)
