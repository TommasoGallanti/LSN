CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -L/opt/homebrew/opt/armadillo/lib -larmadillo
all: main.exe
main.exe : position.o main.o random.o population.o
	$(CC) position.o random.o population.o main.o -o main.exe $(AFLAGS)
main.o : main.cpp position.h population.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
position.o : position.cpp position.h random.h 
	$(CC) -c position.cpp -o position.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
population.o : population.cpp population.h random.h
	$(CC) -c population.cpp -o population.o $(CFLAGS)	
clean :
	rm *.o main.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*