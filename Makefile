FC = gfortran
FFLAGS = -O3 -Wall -Wextra -std=f2018 -fopenmp
SRC_DIR = src
BIN_DIR = bin

# Source files
SRCS = $(SRC_DIR)/kinds.f90 \
       $(SRC_DIR)/particle.f90 \
       $(SRC_DIR)/tree.f90 \
       $(SRC_DIR)/physics.f90 \
       $(SRC_DIR)/collisions.f90 \
       $(SRC_DIR)/initial_conditions.f90 \
       $(SRC_DIR)/main.f90

# Object files (in order of dependency)
OBJS = kinds.o particle.o tree.o physics.o collisions.o initial_conditions.o main.o

TARGET = $(BIN_DIR)/celestia

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) -o $@ $^

# Dependencies
kinds.o: $(SRC_DIR)/kinds.f90
	$(FC) $(FFLAGS) -c $<

particle.o: $(SRC_DIR)/particle.f90 kinds.o
	$(FC) $(FFLAGS) -c $<

tree.o: $(SRC_DIR)/tree.f90 kinds.o particle.o
	$(FC) $(FFLAGS) -c $<

physics.o: $(SRC_DIR)/physics.f90 kinds.o particle.o tree.o
	$(FC) $(FFLAGS) -c $<

collisions.o: $(SRC_DIR)/collisions.f90 kinds.o particle.o
	$(FC) $(FFLAGS) -c $<

initial_conditions.o: $(SRC_DIR)/initial_conditions.f90 kinds.o particle.o
	$(FC) $(FFLAGS) -c $<

main.o: $(SRC_DIR)/main.f90 kinds.o particle.o physics.o collisions.o initial_conditions.o
	$(FC) $(FFLAGS) -c $<

run: all
	./$(TARGET)

clean:
	rm -f *.o *.mod $(TARGET)
