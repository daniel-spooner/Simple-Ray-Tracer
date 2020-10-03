
all: rayT

rayT:	main.o ray.o primitives.o light.o
	g++ main.o ray.o primitives.o light.o -lglut -lGL -lGLU -lGLEW -o rayT

main.o: main.cpp ray.h primitives.h light.h
	g++ -c main.cpp	

ray.o: ray.cpp ray.h
	g++ -c ray.cpp

primitives.o: primitives.cpp primitives.h
	g++ -c primitives.cpp

light.o: light.cpp light.h
	g++ -c light.cpp