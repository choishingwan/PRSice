PRSice_debug: src/*.c*
	g++ -std=c++11 -g -I inc/ -isystem lib/  -pthread src/*.c* -o PRSice_debug
