.PHONY : all clean run
all : run post

run : 
	mpic++ -std=c++11 -fpermissive simulation.cpp -o run.exe
	./run.exe

post:
	python plot.py

clean : 
	rm -rf  *.dat *.pdf run.exe