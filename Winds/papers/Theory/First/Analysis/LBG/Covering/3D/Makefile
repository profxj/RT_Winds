#- 
#- File    : Makefile
#- ------------------
#- Created : Mon Feb 11 17:33:35 2008
#- 
#- $Header: /home/jonathan/xavier/X_CVS/radtrans/Winds/papers/Theory/First/Analysis/LBG/Covering/3D/Makefile,v 1.2 2010-08-20 21:00:54 xavier Exp $
#- 

TARGET = lbg_covering_3D

OBJECTS =  Lua.o main_lbg.o transport_lbg.o anvoigt.o locate_array.o cdf.o spectrum.o model_lbg.o

# additional stuff for pleiades
#LUA_BASE  = /ibrixfs/kasen/ddmc/ddmc_sn/lua-5.1.3
#LUA_SRC   = $(LUA_BASE)/src
#LIB_LUA_A = $(LUA_SRC)/liblua.a
#GSLI  = /home/kasen/local/include
#GSLL  = /home/kasen/local/lib


# additional stuff for macs
LUA_BASE  = /usr/local/
LUA_SRC   = $(LUA_BASE)/lib
LIB_LUA_A = $(LUA_SRC)/liblua.a
GSLI  = /usr/local/include
GSLL  = /usr/local/lib


CXX   = mpicxx -g  -DMPICH_IGNORE_CXX_SEEK
CC    = mpicxx
CCOPT = -I$(GSLI) -I$(LUA_SRC)
CLOPT = -L$(GSLL) -L$(LUA_SRC) -llua -lgsl -lgslcblas


$(TARGET): $(OBJECTS)
	$(CXX) -o $(TARGET) $(OBJECTS) $(CLOPT)

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
