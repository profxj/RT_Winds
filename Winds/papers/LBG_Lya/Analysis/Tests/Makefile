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



TARGET = test_static_sphere
TARGET2 = test_static_sphere2
#TARGET2 = outflow_fiducial_feii

OBJECTS = test_static_sphere.o miki_voigt.o locate_array.o cdf.o spectrum.o 
OBJECTS2 = test_static_sphere2.o miki_voigt.o locate_array.o cdf.o spectrum.o 


all: test_static_sphere test_static_sphere2

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

test_static_sphere: $(OBJECTS) 
	$(CXX) -o $(TARGET) $(OBJECTS) $(CLOPT)

test_static_sphere2: $(OBJECTS2) 
	$(CXX) -o $(TARGET2) $(OBJECTS2) $(CLOPT)

