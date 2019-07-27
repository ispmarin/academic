CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		src/smeagol.o

LIBS =

TARGET =	src/smeagol

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
