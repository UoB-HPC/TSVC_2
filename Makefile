
COMPILER=GNU

ifndef COMPILER
$(error No compiler specified!)
endif

ifndef SRC_DIR
$(warning No SRC_DIR is specified, using current dir as src dir!)
SRC_DIR=./src
endif

ifndef KERNEL_NAME
$(warning No kernel name defined, building for all!)
KERNEL_NAME=ALL_KERNEL
endif

all:
	@$(MAKE) -C $(SRC_DIR) COMPILER=$(COMPILER)                NO_OMP=1 SUFFIX=_default KERNEL_NAME=$(KERNEL_NAME)
	@$(MAKE) -C $(SRC_DIR) COMPILER=$(COMPILER) FAST_MATH=1    NO_OMP=1 SUFFIX=_relaxed KERNEL_NAME=$(KERNEL_NAME)
	@$(MAKE) -C $(SRC_DIR) COMPILER=$(COMPILER) PRECISE_MATH=1 NO_OMP=1 SUFFIX=_precise KERNEL_NAME=$(KERNEL_NAME)

.PHONY: clean
clean:
	$(MAKE) -C $(SRC_DIR) COMPILER=$(COMPILER) KERNEL_NAME=$(KERNEL_NAME) clean

