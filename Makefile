default: cLJA

cLJA: consensus_asm.cpp dot_graph.o graph_simplification.o
	g++ -O3 -g -o cLJA consensus_asm.cpp edlib/src/edlib.cpp dot_graph.o graph_simplification.o -lpthread -fopenmp -fconcepts -I edlib/include

dot_graph.o: dot_graph.cpp dot_graph.hpp
	g++ -O3 -g -c dot_graph.cpp edlib/src/edlib.cpp -fopenmp -fconcepts -I edlib/include
	
graph_simplification.o: graph_simplification.cpp
	g++ -O3 -g -c graph_simplification.cpp edlib/src/edlib.cpp -I edlib/include -fconcepts

clean:
	rm -f cLJA dot_graph.o graph_simplification.o