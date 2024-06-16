CC = g++
CFLAGS = -Wall -O3 --std=c++11


EXECUTABLES = main2.1.exe main2.2.exe

all: 
		make $(EXECUTABLES)



%.exe : %.o random.o
	$(CC) $^ -o $@

%.o : %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)


random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main2.1.exe main2.2.exe