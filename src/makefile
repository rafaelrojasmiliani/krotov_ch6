
SRCGSLEVENT= gslodeeventstep.cpp gslodeeventrootridder.cpp \
			 gslodeeventhermite.cpp gslodeeventrootbrent.cpp indexvector.cpp

SRCSOLVER = solver.cpp solveraux.cpp

SRCKRVITER= iterator.cpp  iteratordynsys.cpp

SRCKRVAFFN= affinesys.cpp affinesysproy.cpp affinesyscontrol.cpp \
			affinesysaux.cpp

SRCSOLCONT=  csolcont.cpp  csolcont.h  csolcontio.cpp csolcontiofile.cpp

Warning=-Wuninitialized  -Wcast-qual -Wsuggest-attribute=const  -std=c++11 -DNDEBUG


all: 	libiteratorgdb.a libgslodeeventgdb.a \
		libaffinesysgdb.a libsolcontgdb.a libsolvergdb.a \
		libiterator.a libgslodeevent.a \
		libaffinesys.a libsolcont.a libsolver.a

	ar cr libkrotovch6gdb.a 	$(SRCKRVAFFN:%.cpp=%gdb.o) \
								$(SRCSOLVER:%.cpp=%gdb.o) \
								$(SRCKRVITER:%.cpp=%gdb.o) \
								$(SRCGSLEVENT:%.cpp=%gdb.o) \
								$(SRCSOLCONT:%.cpp=%gdb.o) 

	ar cr libkrotovch6.a 	$(SRCKRVAFFN:%.cpp=%.o) \
							$(SRCSOLVER:%.cpp=%.o) \
							$(SRCKRVITER:%.cpp=%.o) \
							$(SRCGSLEVENT:%.cpp=%.o) \
							$(SRCSOLCONT:%.cpp=%.o) 

	cp iterator.h 			~/programacion/include/	
	cp solver.h 			~/programacion/include/	
	cp affinesys.h 			~/programacion/include/	
	cp gslodeevent.h 		~/programacion/include/	
	cp csolcont.h 			~/programacion/include/	
	cp indexvector.h 		~/programacion/include/	
	mv libkrotovch6gdb.a 	~/programacion/lib
	mv libkrotovch6.a 		~/programacion/lib
	

opt: 	libiterator.a libgslodeevent.a \
		libaffinesys.a libsolcont.a libsolver.a
	g++ $(Warning)  -O3 test.cpp -L. \
			 -laffinesys -lgslodeevent \
				-lsolver -literator   -lsolcont  -lgsl -lblas -o test

ld: 	libiteratorgdb.a libgslodeeventgdb.a \
		libaffinesysgdb.a libsolcontgdb.a libsolvergdb.a
	g++ $(Warning)  -g3 testloaddata.cpp -L. \
					 -laffinesysgdb \
					-lgslodeeventgdb -lsolvergdb -literatorgdb -lsolcontgdb \
					-lgsl -lblas -o test
r: 	libiteratorgdb.a libgslodeeventgdb.a \
		libaffinesysgdb.a libsolcontgdb.a libsolvergdb.a
	g++ $(Warning)  -g3 testr.cpp -L. \
					 -laffinesysgdb \
					-lgslodeeventgdb -lsolvergdb -literatorgdb -lsolcontgdb \
					-lgsl -lblas -o test

optld: 	libiterator.a libgslodeevent.a \
		libaffinesys.a libsolcont.a libsolver.a
	g++ $(Warning)  -O3 testloaddata.cpp -L. \
			 -laffinesys -lgslodeevent \
				-lsolver -literator   -lsolcont  -lgsl -lblas -o test

optr: 	libiterator.a libgslodeevent.a \
		libaffinesys.a libsolcont.a libsolver.a
	g++ $(Warning)  -O3 testr.cpp -L. \
			 -laffinesys -lgslodeevent \
				-lsolver -literator   -lsolcont  -lgsl -lblas -o test
optt: 	libiterator.a libgslodeevent.a \
		libaffinesys.a libsolcont.a libsolver.a
	g++ $(Warning)  -O3 testthread.cpp -L. \
			 -laffinesys -lgslodeevent \
				-lsolver -literator   -lsolcont  -lgsl -lblas  -pthread -o test

thread: 	libiteratorgdb.a libgslodeeventgdb.a \
		libaffinesysgdb.a libsolcontgdb.a libsolvergdb.a
	g++ $(Warning)  -g3 testthread.cpp -L. \
			 -laffinesysgdb -lgslodeeventgdb \
				-lsolvergdb -literatorgdb   -lsolcontgdb  -lgsl -lblas  -pthread -o test
#{{{ solcont gdb
libsolcontgdb.a: $(SRCSOLCONT:%.cpp=%gdb.o) 
	ar cr $@ $(SRCSOLCONT:%.cpp=%gdb.o)
csolcontgdb.o: csolcont.cpp csolcont.h
	g++ $(Warning) -g3 $< -c -o $@
csolcont%gdb.o: csolcont%.cpp csolcont.h
	g++ $(Warning) -g3 $< -c -o $@
#}}}

#{{{ solver gdb
libsolvergdb.a: $(SRCSOLVER:%.cpp=%gdb.o) 
	ar cr $@ $(SRCSOLVER:%.cpp=%gdb.o) 
solvergdb.o: solver.cpp solver.h iterator.h
	g++ $(Warning) -g3 $< -c -o $@
solver%gdb.o: solver%.cpp solver.h iterator.h
	g++ $(Warning) -g3 $< -c -o $@
#}}}

#{{{ iteratorgdb
libiteratorgdb.a: $(SRCKRVITER:%.cpp=%gdb.o) 
	ar cr $@ $(SRCKRVITER:%.cpp=%gdb.o)

iteratorgdb.o: iterator.cpp iterator.h libsolcontgdb.a
	g++ $(Warning) -g3 $< -c -o $@
iterator%gdb.o: iterator%.cpp iterator.h libsolcontgdb.a
	g++ $(Warning) -g3 $< -c -o $@
#}}}

#{{{ continuousrkgdb
libgslodeeventgdb.a: $(SRCGSLEVENT:%.cpp=%gdb.o) 
	ar cr $@  $(SRCGSLEVENT:%.cpp=%gdb.o)
indexvectorgdb.o: indexvector.cpp indexvector.h
	g++ $(Warning) -std=c++11 -g3 $< -c -o $@
gslodeevent%gdb.o: gslodeevent%.cpp gslodeevent.h
	g++ $(Warning) -std=c++11 -g3 $< -c -o $@
#}}}


#{{{ affinesysgdb 
libaffinesysgdb.a: $(SRCKRVAFFN:%.cpp=%gdb.o)
	ar cr $@  $(SRCKRVAFFN:%.cpp=%gdb.o)
affinesysgdb.o: affinesys.cpp affinesys.h solver.h iterator.h
	g++ $(Warning) -g3 $< -c -o $@
affinesys%gdb.o: affinesys%.cpp affinesys.h solver.h iterator.h
	g++ $(Warning) -g3 $< -c -o $@
#}}}

##------------------------------- opt ------

#{{{ solcont opt
libsolcont.a: $(SRCSOLCONT:%.cpp=%.o) 
	ar cr $@ $(SRCSOLCONT:%.cpp=%.o)
csolcont.o: csolcont.cpp csolcont.h
	g++ $(Warning) -O3 $< -c -o $@
csolcont%.o: csolcont%.cpp csolcont.h
	g++ $(Warning) -O3 $< -c -o $@
#}}}

#{{{ solver opt
libsolver.a: $(SRCSOLVER:%.cpp=%.o) 
	ar cr $@ $(SRCSOLVER:%.cpp=%.o)
solver.o: solver.cpp solver.h
	g++ $(Warning) -O3 $< -c -o $@
solver%.o: solver%.cpp solver.h
	g++ $(Warning) -O3 $< -c -o $@
#}}}
#{{{ iterator opt
libiterator.a: $(SRCKRVITER:%.cpp=%.o) 
	ar cr $@ $(SRCKRVITER:%.cpp=%.o)

iterator.o: iterator.cpp iterator.h
	g++ $(Warning) -O3 $< -c -o $@
iterator%.o: iterator%.cpp iterator.h
	g++ $(Warning) -O3 $< -c -o $@
#}}}

#{{{ continuousrk opt
libgslodeevent.a: $(SRCGSLEVENT:%.cpp=%.o) 
	ar cr $@  $(SRCGSLEVENT:%.cpp=%.o)

indexvector.o: indexvector.cpp indexvector.h
	g++ $(Warning) -std=c++11 -O3 $< -c -o $@
gslodeevent%.o: gslodeevent%.cpp gslodeevent.h
	g++ $(Warning) -std=c++11 -O3 $< -c -o $@
#}}}

#{{{ affinesys opt
libaffinesys.a: $(SRCKRVAFFN:%.cpp=%.o)
	ar cr $@  $(SRCKRVAFFN:%.cpp=%.o)

affinesys.o: affinesys.cpp affinesys.h solver.h
	g++ $(Warning) -O3 $< -c -o $@
affinesys%.o: affinesys%.cpp affinesys.h solver.h
	g++ $(Warning) -O3 $< -c -o $@
#}}}

clean:
	-rm *.o
	-rm lib*.a


cleanall:
	-rm *.o
	-rm lib*.a
	-rm *.dat


doc:
	doxygen  && cd latex && make 

#
#iterator.o: $(SRCKRVITER) 
#	g++ $(Warning) -O3 $< -c $@
#
#
#contrkgdb.o: $(SRCGSLEVENT)
#	g++ $(Warning) -g3 $< -c $@
#contrk.o: $(SRCGSLEVENT)
#	g++ $(Warning) -O3 $< -c $@
#
#affinesysgdb.o: $(SRCKRVAFFN)
#	g++ $(Warning) -g3 $< -c $@
#affinesys.o: $(SRCKRVAFFN)
#	g++ $(Warning) -O3 $< -c $@
#
#solver:
#	g++ $(Warning)  -O3 testsolver.cpp iterator.cpp iteratorio.cpp iteratordynsys.cpp  -lgsl -lblas -o test
#
#ehmb:
#	g++ $(Warning)  -g3 testehmb.cpp iterator.cpp iteratorio.cpp iteratordynsys.cpp  affinesysstep.cpp affinesyssteprk2.cpp -lgsl -lblas -o test
#
