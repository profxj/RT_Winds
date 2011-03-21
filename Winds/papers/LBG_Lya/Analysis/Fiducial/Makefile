#- 
#- File    : Makefile
#- ------------------
#- Created : Mon Feb 11 17:33:35 2008
#- 
#- $Header: /home/jonathan/xavier/X_CVS/radtrans/Winds/papers/Theory/First/Analysis/Fiducial/Makefile,v 1.2 2010-06-11 13:14:38 xavier Exp $
#- 

# For Dan's Mac
 CXX = mpicxx
 CCOPT = $(GSL) 
 CLOPT = $(GSL)  -lgsl 



TARGET = outflow_fiducial_lya
#TARGET2 = outflow_fiducial_feii

OBJECTS = outflow_fiducial_lya.o jxp_voigt.o locate_array.o cdf.o spectrum.o 
#OBJECTS2 = outflow_fiducial_feii.o voigt.o locate_array.o cdf.o spectrum.o 


all: lya

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

lya: $(OBJECTS) 
	$(CXX) -o $(TARGET) $(OBJECTS) $(CLOPT)

