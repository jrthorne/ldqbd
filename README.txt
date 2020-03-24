1. Compile modules first
$ gcc -c matutil.c -o matutil.o
$ gcc -c ldqbd.c -o ldqbd.o
$ gcc -c qgen.c -o qgen.o

2. Then compile the main code
$ gcc matmain.c ldqbd.o matutil.o qgen.o -o matmain

3. Then run with
$ ./matmain

NOTE: I based my honours thesis pretty heavily on Les Bright's PhD Thesis, so page numbers referred to are generally to do with his thesis. I still have a hard copy of this (I think), but Peter or  will probably still have a hard copy too.