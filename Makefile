PROGRAM       = BDT

LIBDIRARCH	  = lib

CXX           = /usr/bin/g++
CXXFLAGS      = -Wall -fPIC -g

LD            = g++
SOFLAGS       = -shared
CXXFLAGS     += $(shell root-config --cflags)
LIBS          = $(shell root-config --libs)

#LIBS		 += -L$(PYTHIA8)/$(LIBDIRARCH) -lpythia8
#INCS		  = -I$(PYTHIA8)/include

#INCS += $(shell $(FASTJET)/bin/fastjet-config --cxxflags )         # /fastjet
#LIBS += $(shell $(FASTJET)/bin/fastjet-config --libs --plugins )   # $(INCS)

HDRSDICT 	= 

HDRS		+= $(HDRSDICT) 
SRCS 		= $(HDRS:.h=.C)
OBJS		= $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
		$(CXX) $(CXXFLAGS) $(OBJS) $(PROGRAM).C $(LIBS) $(INCS) \
		-o $(PROGRAM)
		@echo "$(PROGRAM) done"
$(OBJS): %.cc %.h
		$(CXX) $(CXXFLAGS) $(LIBS) $(INCS) -o $@ $<

clean:
		rm -f $(OBJS) core *Dict*  $(PROGRAM).o $(PROGRAM)

headerDict.C:	$(HDRSDICT)
		@echo "Generating dictionary ..."
		@rm -f headerDict.C headerDict.h
		@rootcint headerDict.C -c $(HDRSDICT)

