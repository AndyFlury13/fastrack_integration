CC=g++
CFLAGS= -O3 -std=c++11 -fopenmp -mfma -msse2 -c -Wall -Wextra -fPIC #-ffast-math

LDFLAGS += -lgomp -lrt

SOURCES_FFF = AlgorithmParameters.cpp Geometry.cpp LayerLinker.cpp FieldModel.cpp Track.cpp
SOURCES_FFF += HitBin.cpp HitLayer.cpp HitBank.cpp ClusterMaker.cpp
SOURCES_FFF += TrackingFilter.cpp FastTrackExtrapolator.cpp FastTrackFitter.cpp
SOURCES_FFF += TrackFinder.cpp ModelClass.cpp 

SOURCE_FFF = Model.cpp
OBJECTS_FFF=$(SOURCES_FFF:.cpp=.o)
OBJECT_FFF=$(SOURCE_FFF:.cpp=.o)
LIBRARY_FFF=libModel.so

pro: $(SOURCES_FFF) $(SOURCE_FFF) $(LIBRARY_FFF)

$(LIBRARY_FFF): $(OBJECTS_FFF)
	$(CC) $(CFLAGS) -c $(SOURCE_FFF) -o $(OBJECT_FFF)
	$(CC) $(LDFLAGS) -shared -Wl,-soname,$(LIBRARY_FFF) -o $(LIBRARY_FFF) $(OBJECT_FFF) $(OBJECTS_FFF)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@


clean:
	rm *.o *.so
