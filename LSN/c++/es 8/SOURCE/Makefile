CC = g++
CFLAGS = -O3 --std=c++11
AFLAGS = -L/opt/homebrew/opt/armadillo/lib -larmadillo
all: main.exe
main.exe : varMC.o main.o random.o
	$(CC) varMC.o random.o main.o -o main.exe $(AFLAGS)
main.o : main.cpp varMC.h 
	$(CC) -c main.cpp -o main.o $(CFLAGS)
varMC.o : varMC.cpp varMC.h random.h 
	$(CC) -c varMC.cpp -o varMC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe
remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*