TARGET = test.exe
OBJS = main.o
CC = g++

# boost include for mac
INC=-I/opt/homebrew/Cellar/boost/1.80.0/include

CFLAGS = -c -Wall -g -std=c++11
LFLAGS = -Wall -g
#CFLAGS = -c -Wall -O3 -DNDEBUG -std=c++11
#LFLAGS = -Wall -O3 -DNDEBUG

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(TARGET)

main.o: main.cpp vec3.h config_file.h systemEDBD.h initialize_positions.h
	$(CC) $(CFLAGS) $(INC) main.cpp


.PHONY: clean
clean:
	rm -f  $(OBJS) $(TARGET) 

.PHONY: cleanObject
cleanObject:
	rm -f  $(OBJS)

