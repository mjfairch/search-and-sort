all:
	clang -Wall -O3 SearchAndSort.c -o sort

test:
	./sort
