all: decompose

decompose: decompose.cc
	g++ decompose.cc -o decompose

test: decompose
	./decompose > test-output.txt
	diff test-output.txt known-good.txt
