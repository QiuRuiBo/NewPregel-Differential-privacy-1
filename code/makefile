DEPS=src/GraphEnginePregelAggregateDC.cpp GeoDistrWorkflow.cpp src/GraphConstruction.cpp src/Utils.cpp src/countTriangles/Graph1.cpp src/countTriangles/Subgraph3Counter.cpp src/Count_triangles.cpp FindVertexRandomNeighbor.cpp src/Check_connectivity/Count_subgraphs.cpp

OSFLAG = $(shell uname -s)
ifeq ($(OSFLAG),Darwin)
	CC =clang++
	CFLAGS =--std=c++11 -Xpreprocessor -fopenmp -lomp -g
else ifeq ($(OSFLAG),Linux)
	CC=/usr/bin/g++-7
	CFLAGS=-I ../boost_1_50_0 -std=c++11  -fopenmp -g
endif
PrivacyRun: $(DEPS)
	$(CC) $(CFLAGS) $(DEPS) -o PrivacyRun
#	$(CC) $(CFLAGS) $(DEPS) -o PrivacyRun
