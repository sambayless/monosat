################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include src/monosat/utils/subdir.mk
-include src/monosat/simp/subdir.mk
-include src/monosat/pb/ADTs/subdir.mk
-include src/monosat/pb/subdir.mk
-include src/monosat/graph/subdir.mk
-include src/monosat/geometry/polypartition/subdir.mk
-include src/monosat/geometry/cevans/subdir.mk
-include src/monosat/geometry/bitmap/subdir.mk
-include src/monosat/geometry/subdir.mk
-include src/monosat/fsm/subdir.mk
-include src/monosat/dgl/alg/subdir.mk
-include src/monosat/dgl/subdir.mk
-include src/monosat/core/subdir.mk
-include src/monosat/api/subdir.mk
-include src/monosat/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: libmonosat_dbg.so

# Tool invocations
libmonosat_dbg.so: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++  -shared -o "libmonosat_dbg.so" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(LIBRARIES)$(CC_DEPS)$(C++_DEPS)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS) libmonosat_dbg.so
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
