CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -L/opt/homebrew/opt/armadillo/lib -larmadillo
all: simulator.exe
simulator.exe : system.o main.o particle.o random.o
	$(CC) system.o particle.o random.o main.o -o simulator.exe $(AFLAGS)
main.o : main.cpp system.h particle.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
system.o : system.cpp system.h random.h particle.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)
particle.o : particle.cpp particle.h
	$(CC) -c particle.cpp -o particle.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o simulator.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
