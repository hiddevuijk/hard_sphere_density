TARGET = test.exe
OBJS = main.o xyz.o
CC = g++
#CFLAGS = -c -Wall -g -std=c++11
#LFLAGS = -Wall -g
CFLAGS = -c -Wall -O3 -DNDEBUG -std=c++11
LFLAGS = -Wall -O3 -DNDEBUG

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS)  $(OBJS) -o $(TARGET)

xyz.o: xyz.cpp xyz.h
	$(CC) $(CFLAGS) xyz.cpp

main.o: main.cpp xyz.h
	$(CC) $(CFLAGS) main.cpp


.PHONY: clean
clean:
	rm -f  $(OBJS) $(TARGET) 

.PHONY: cleanObject
cleanObject:
	rm -f  $(OBJS)

