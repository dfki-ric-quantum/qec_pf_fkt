SHELL     = /bin/bash
CXX       = g++
CXXFLAGS  = -O3 -Wall -W -pedantic
LIBS      = -lgslcblas -lgsl

BUILD_DIR = ../../build/generator_random_bond
PROGNAME  = $(BUILD_DIR)/isingGeneratorRandomBond

SRCS      = main.cc
OBJS      = $(SRCS:%.cc=$(BUILD_DIR)/%.o)

all: $(PROGNAME)

# Ensure build directory exists
$(BUILD_DIR):
	mkdir -p $@

# Compile object files in build dir
$(BUILD_DIR)/%.o: %.cc | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link final binary into build dir
$(PROGNAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	@rm -f $(BUILD_DIR)/*.o $(PROGNAME)

.PHONY: all clean
