cd Documents/Github/SimplexSIS/src
g++ -c SimplexContagion.cpp -llapack -lblas
g++ -shared -o SimplexContagion.dll SimplexContagion.o
