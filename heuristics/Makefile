TWOUP    = ..
CPP      = g++
CARGS    = -m64 -g -Wall
HEADERS =  my_struct.h load_input.h priority_queue.h preprocessing_by_collapse.h preprocessing_by_components.h propagation.h greedy_randomized.h weighted_greedy_randomized.h random_plus_greedy.h sampled_greedy.h local_search.h grasp.h
SOURCES =  my_struct.cpp load_input.cpp priority_queue.cpp preprocessing_by_collapse.cpp preprocessing_by_components.cpp propagation.cpp greedy_randomized.cpp weighted_greedy_randomized.cpp random_plus_greedy.cpp sampled_greedy.cpp local_search.cpp grasp.cpp grasp_pa_main.cpp

grasp_pa: $(HEADERS) $(SOURCES)
	$(CPP) $(CARGS) -o  $@ $(SOURCES)

run_grasp_pa: grasp_pa
	./pli_pa $(args)

clean:
	rm -rf *.o