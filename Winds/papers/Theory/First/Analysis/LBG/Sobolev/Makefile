#- 
#- File    : Makefile
#- ------------------
#- Created : Mon Feb 11 17:33:35 2008
#- 
#- $Header: /home/jonathan/xavier/X_CVS/radtrans/Winds/papers/Theory/First/Analysis/LBG/Sobolev/Makefile,v 1.1 2010-06-24 12:27:49 xavier Exp $
#- 

# For Dan's Mac
 CXX = mpicxx
 CCOPT = $(GSL) 
 CLOPT = $(GSL)  -lgsl 



TARGET = outflow_lbg_sobolev_mgii
TARGET2 = outflow_lbg_sobolev_feii

OBJECTS = outflow_lbg_sobolev_mgii.o voigt.o locate_array.o cdf.o spectrum.o 
OBJECTS2 = outflow_lbg_sobolev_feii.o voigt.o locate_array.o cdf.o spectrum.o 


all: mgii

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

mgii: $(OBJECTS) 
	$(CXX) -o $(TARGET) $(OBJECTS) $(CLOPT)

feii: $(OBJECTS2) 
	$(CXX) -o $(TARGET2) $(OBJECTS2) $(CLOPT)

