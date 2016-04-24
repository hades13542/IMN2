all:
	g++ -o main main.cpp
	./main
	gnuplot gnuplot
git:
	git add --all
	git commit -m "auto-generate commit by makefile"
	git push