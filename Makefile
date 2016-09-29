build : main.o
	gcc main.o -o quadtree -lm
main.o : main.c
	gcc -g -c main.c
clean :
	rm quadtree main.o
