CFLAGS=-Wall -Wextra -O3 #-Werror   
LFLAGS=-lstdc++ -fopenmp -lm -lfftw3 -lfftw3f #-lpng -ljpeg -ltiff 

#object files
OBJ_ICA= bicubic_interpolation.o file.o inverse_compositional_algorithm.o mask.o matrix.o transformation.o zoom.o

OBJ_ESTADEO= color_bicubic_interpolation.o crop_and_zoom.o estadeo.o gaussian_conv_dct.o main.o motion_smoothing.o online_motion_smoothing.o sutherland_hodgman.o utils.o

OBJ= $(OBJ_ICA) $(OBJ_ESTADEO)

#executable files
all: bin obj bin/estadeo bin/generate_graphics

bin:
	mkdir -p bin

obj:
	mkdir -p obj

#generate executables
bin/estadeo: $(addprefix obj/,$(OBJ)) 
	g++ $^ -o $@ $(CFLAGS) $(LFLAGS)

bin/generate_graphics: obj/generate_graphics.o obj/gaussian_conv_dct.o  
	g++ $^ -o $@ $(CFLAGS) $(LFLAGS) 


#compile ica
obj/%.o: src/ica/%.cpp
	g++ -c $< -o $@ $(INCLUDE) $(CFLAGS) $(LFLAGS)

#compile estadeo 
obj/%.o: src/%.cpp
	g++ -c $< -o $@ $(INCLUDE) $(CFLAGS) $(LFLAGS)

clean: 
	rm -f bin/estadeo bin/generate_graphics
	rm -R obj

