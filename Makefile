# Makefile

my_program: main.o main.h
	g++ -std=c++17 main.o -o my_program
    
main.o: main.cpp main.h
	g++ -std=c++17 -c -Wall -Wextra main.cpp 
    
clean:
	rm -f *.o my_program