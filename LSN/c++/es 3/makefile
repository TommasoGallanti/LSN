CC = g++
CFLAGS = -Wall -O3 --std=c++11


EXECUTABLES = main.exe 

all: 
		make $(EXECUTABLES)



%.exe : %.o random.o
	$(CC) $^ -o $@

%.o : %.cpp
	$(CC) -c $< -o $@ $(CFLAGS)


random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main.exe