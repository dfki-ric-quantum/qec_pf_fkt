SHELL      = /bin/bash
CXX        = g++
CXXFLAGS   = -m64 -O3 -Wall -W -pedantic
LIBS       = -lgslcblas -lgsl -lgmp -lgmpxx

BUILD_DIR  = ../../build/Z_to_txt
PROGNAME   = $(BUILD_DIR)/isingZToTxt

SRCS       = main.cc FINDmatrix.cc Sample.cc exp_log.cc
OBJS       = $(SRCS:%.cc=$(BUILD_DIR)/%.o)

all: $(PROGNAME)

# Create build directory if needed
$(BUILD_DIR):
	@mkdir -p $@

# Build the binary inside build dir
$(PROGNAME): $(BUILD_DIR) $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

# Compile source files to object files inside build dir
$(BUILD_DIR)/%.o: %.cc | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@rm -f $(BUILD_DIR)/*.o $(PROGNAME)

.PHONY: all clean
