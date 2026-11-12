
include ../makefiles/Makefile.defs
# Defines the list of tests
include loops.mk

COMMON_FILES=dummy.o common.o

BIN_DIR=../bin/$(COMPILER)

libs=-lm

all : $(BIN_DIR)/tsvc_vec $(BIN_DIR)/tsvc_novec

$(BIN_DIR)/tsvc_% : tsvc_%.o $(foreach t,$(TESTS), loops/$(t)_%.o) $(COMMON_FILES)
	@mkdir -p $(BIN_DIR)
	$(CC) $^ $(libs) -o $@$(SUFFIX)

%_novec.o:%.c .FORCE
	$(CC) -I . $(flags) $(novecflags) -c -o $@ $<

%_vec.o:%.c .FORCE
	$(CC) -I . $(flags) $(vecflags) -c -o $@ $<

# Compiled with IEEE math flags to make initialisation/checksum calculations more reproducable
# Be wary that NVHPC doesn't respect IEEE properly with vectorisation
%.o:%.c .FORCE
	$(CC) $(helper_flags) $(vecflags) $(ieee_math_flags) -c -o $@ $<

.PHONY: .FORCE clean
.FORCE:
	@true

clean :
	$(RM) -f *.o src/*.o *.s

