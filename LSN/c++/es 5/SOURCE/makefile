CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -L/opt/homebrew/opt/armadillo/lib -larmadillo
all: main.exe
main.exe : main.o particle.o metropolis.o random.o 
	$(CC) particle.o metropolis.o random.o main.o -o main.exe $(AFLAGS)
main.o : main.cpp particle.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
metropolis.o : metropolis.cpp metropolis.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
particle.o : particle.cpp particle.h
	$(CC) -c particle.cpp -o particle.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*