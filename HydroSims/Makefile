#- 
#- File    : Makefile
#- ------------------
#- Created : Mon Feb 11 17:33:35 2008
#- 
#- $Header: /home/jonathan/xavier/X_CVS/radtrans/Winds/Grid/Transfer/Makefile,v 1.1 2010-07-30 16:45:53 xavier Exp $
#- 

TARGET = a.out

OBJECTS =  Lua.o main.o transport.o anvoigt.o locate_array.o cdf.o spectrum.o model.o uv_transport.o

# additional stuff for pleiades
#LUA_BASE  = /ibrixfs/kasen/ddmc/ddmc_sn/lua-5.1.3
#LUA_SRC   = $(LUA_BASE)/src
LUA_BASE  = /usr/local/lib
LUA_SRC   = $(LUA_BASE)
LIB_LUA_A = $(LUA_SRC)/liblua.a
GSLI  = /usr/local/include/gsl
GSLL  = /usr/local/lib
GDL = /usr/lib/x86_64-linux-gnu


# additional stuff for macs
#LUA_BASE  = /usr/local/
#LUA_SRC   = $(LUA_BASE)/lib
#LIB_LUA_A = $(LUA_SRC)/liblua.a
#GSLI  = /usr/local/include
#GSLL  = /usr/local/lib


CXX   = mpicxx -g  -DMPICH_IGNORE_CXX_SEEK
CC    = mpicxx
CCOPT = -I$(GSLI) -I$(LUA_SRC)
CLOPT = -L$(GSLL) -L$(LUA_SRC) -llua -lgsl -lgslcblas
#LDLOPT = -L$(GDL) -ldl
LDLOPT = /usr/lib/x86_64-linux-gnu/libdl.so


$(TARGET): $(OBJECTS)
	$(CXX) -o $(TARGET) $(OBJECTS) $(CLOPT) $(LDLOPT)

.cc.o: 
	$(CXX) $(CCOPT) -c -o $@ $<

.c.o: 
	cc -I. -c -o $@ $<

clean:
	rm -f $(OBJECTS)

logclean:
	rm -f output.*.log

realclean: clean
	rm -f core $(TARGET)
