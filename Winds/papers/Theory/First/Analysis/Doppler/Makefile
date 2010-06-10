#- 
#- File    : Makefile
#- ------------------
#- Created : Mon Feb 11 17:33:35 2008
#- 
#- $Header: /home/jonathan/xavier/X_CVS/radtrans/Winds/papers/Theory/First/Analysis/Doppler/Makefile,v 1.1 2010-06-10 19:17:39 xavier Exp $
#- 

# For Dan's Mac
 CXX = mpicxx
 CCOPT = $(GSL) 
 CLOPT = $(GSL)  -lgsl 


TARGET = outflow_doppler_mgii

OBJECTS = outflow_doppler_mgii.o voigt.o locate_array.o cdf.o spectrum.o 


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
