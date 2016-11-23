PRSice:
	g++ -std=c++11 -I inc/ -isystem lib/  -DNDEBUG -O2 -lpthread src/*.cpp -o PRSice
