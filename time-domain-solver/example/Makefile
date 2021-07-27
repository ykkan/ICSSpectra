.PHONY : all clean run
all : run post

run : 
	mpic++ -std=c++11 -fpermissive simulation.cpp -o run
	./run

post:
	julia spectralangular_post.jl
	python plot.py

clean : 
	rm -rf  rawdata/ *.dat *.pdf run