CC = g++
AR = ar
CFLAGS   = -O3 -Wall -c -fmessage-length=0 -Wno-unused -Wno-unused-local-typedefs
INCLUDES = /usr/local/boost ./include
SOURCES  = src/GeomException.cpp \
	   src/random/LCRand.cpp \
	   src/random/MTRand.cpp \
	   src/random/Random.cpp \
	   src/random/RandomTools.cpp

OBJECTS  = $(addprefix build/, $(notdir $(SOURCES:.cpp=.o)))
LIB      = lib/libgeomc.a
INCDIR   = /opt/local/include
LIBDIR   = /opt/local/lib

all: lib
	echo done.

lib: $(OBJECTS)
	$(AR) rs $(LIB) $(OBJECTS)

build/GeomException.o : src/GeomException.cpp
	$(CC) $(CFLAGS) $(addprefix -I, $(INCLUDES)) -o build/GeomException.o src/GeomException.cpp

build/%.o : src/random/%.cpp
	$(CC) $(CFLAGS) $(addprefix -I, $(INCLUDES)) -o $@ $<

install:
	mkdir -p $(INCDIR)
	cp -rf ./include/geomc $(INCDIR)
	cp -rf $(LIB) $(LIBDIR)

clean:
	rm -f ./build/*.o
	rm -f ./lib/*.a

uninstall:
	rm -rf $(INCDIR)/geomc
	rm -f $(LIBDIR)/$(LIB)
