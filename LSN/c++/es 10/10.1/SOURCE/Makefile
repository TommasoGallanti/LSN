CXX = mpic++
CFLAGS = -O3 --std=c++11
AFLAGS = -L/opt/homebrew/opt/armadillo/lib -larmadillo
LDFLAGS = -L/opt/homebrew/Cellar/libevent/2.1.12_1/lib
INCFLAGS = -I/opt/homebrew/Cellar/libevent/2.1.12_1/include

all: main.exe

main.exe : position.o main.o random.o population.o
	$(CXX) position.o random.o population.o main.o -o main.exe $(AFLAGS) $(LDFLAGS)

main.o : main.cpp position.h population.h
	$(CXX) -c main.cpp -o main.o $(CFLAGS) $(INCFLAGS)

position.o : position.cpp position.h random.h 
	$(CXX) -c position.cpp -o position.o $(CFLAGS) $(INCFLAGS)

random.o : random.cpp random.h
	$(CXX) -c random.cpp -o random.o $(CFLAGS) $(INCFLAGS)

population.o : population.cpp population.h random.h
	$(CXX) -c population.cpp -o population.o $(CFLAGS) $(INCFLAGS)

clean :
	rm *.o main.exe

remove:
	rm ../OUTPUT/*.*
	rm ../OUTPUT/CONFIG/*.*
