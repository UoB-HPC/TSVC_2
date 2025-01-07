
include ../makefiles/Makefile.defs

FLAGS=

FLAGS+= -D$(KERNEL_NAME)=1

COMMON_FILES=dummy.o common.o

BIN_DIR=../bin/$(COMPILER)/$(KERNEL_NAME)

libs=-lm

MKDIR = mkdir -p

all : $(BIN_DIR)/tsvc_vec $(BIN_DIR)/tsvc_novec

$(BIN_DIR)/tsvc_% : tsvc_%.o $(COMMON_FILES)
	@mkdir -p $(BIN_DIR)
	$(CC) $(addprefix $(BIN_DIR)/, $^) $(libs) $(FLAGS) -o $@$(SUFFIX)

%_novec.o:%.c .FORCE
	@$(MKDIR) $(dir $(BIN_DIR)/$@)
	$(CC) $(flags) $(novecflags) $(FLAGS) -c -o $(BIN_DIR)/$@ $<

%_vec.o:%.c .FORCE
	@$(MKDIR) $(dir $(BIN_DIR)/$@)
	$(CC) $(flags) $(vecflags) $(FLAGS) -c -o $(BIN_DIR)/$@ $<

# Compiled with IEEE math flags to make initialisation/checksum calculations more reproducable
# Be wary that PGI doesn't respect IEEE properly with vectorisation
%.o:%.c .FORCE
	@$(MKDIR) $(dir $(BIN_DIR)/$@)
	$(CC) $(flags) $(vecflags) $(ieee_math_flags) $(FLAGS) -c -o $(BIN_DIR)/$@ $<

.PHONY: .FORCE clean
.FORCE:
	@true

clean :
	rm -f $(BIN_DIR)/*.o $(BIN_DIR)/*.s

