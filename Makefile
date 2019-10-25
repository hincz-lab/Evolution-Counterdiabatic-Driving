# Project-specific settings
PROJECT := n_dimensions
EMP_DIR := ../Empirical/source

# Flags to use regardless of compiler
CFLAGS_all := -Wall -Wno-unused-function -std=c++17 -I$(EMP_DIR)/

# Native compiler information
CXX_nat := g++
CFLAGS_nat := -O3 -DNDEBUG $(CFLAGS_all)
CFLAGS_nat_debug := -g $(CFLAGS_all)

# Emscripten compiler information
CXX_web := emcc
OFLAGS_web_all := -s "EXTRA_EXPORTED_RUNTIME_METHODS=['ccall', 'cwrap', 'writeStringToMemory']" -s TOTAL_MEMORY=67108864 --js-library $(EMP_DIR)/web/library_emp.js --js-library $(EMP_DIR)/web/d3/library_d3.js -s EXPORTED_FUNCTIONS="['_main', '_empCppCallback']" -s DISABLE_EXCEPTION_CATCHING=1 -s NO_EXIT_RUNTIME=1 #--embed-file configs
OFLAGS_web := -Oz -DNDEBUG
OFLAGS_web_debug := -g4 -Oz -pedantic -Wno-dollar-in-identifier-extension

CFLAGS_web := $(CFLAGS_all) $(OFLAGS_web) $(OFLAGS_web_all)
CFLAGS_web_debug := $(CFLAGS_all) $(OFLAGS_web_debug) $(OFLAGS_web_all)


default: nd
native: nd
web: n_dimensions.js
all: 1d nd n_dimensions.js

debug:	CFLAGS_nat := $(CFLAGS_nat_debug)
debug:	nd

# debug-web:	CFLAGS_web := $(CFLAGS_web_debug)
# debug-web:	$(PROJECT).js

# web-debug:	debug-web

nd:	source/$(PROJECT).cc source/$(PROJECT).h
	$(CXX_nat) $(CFLAGS_nat) source/$(PROJECT).cc -o $(PROJECT)
	cp config/NDim.cfg .

1d:	source/ABMtoFP_Evol.c
	$(CXX_nat) source/ABMtoFP_Evol.c -o 1_dimension

n_dimensions.js: source/n_dimensions_web.cc
	$(CXX_web) $(CFLAGS_web) source/n_dimensions_web.cc -o web/n_dimensions.js

clean:
	rm -f $(PROJECT) web/n_dimensions.js web/*.js.map web/*.js.map *~ *.o

# Debugging information
print-%: ; @echo '$(subst ','\'',$*=$($*))'
