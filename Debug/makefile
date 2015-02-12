################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include utils/subdir.mk
-include core/subdir.mk
-include aiger/subdir.mk
-include Static/utils/subdir.mk
-include Static/simp/subdir.mk
-include Static/graph/subdir.mk
-include Static/geometry/polypartition/subdir.mk
-include Static/geometry/cevans/subdir.mk
-include Static/geometry/subdir.mk
-include Static/fsm/subdir.mk
-include Static/dgl/alg/subdir.mk
-include Static/dgl/subdir.mk
-include Static/core/subdir.mk
-include Static/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: modsat_allsat

# Tool invocations
modsat_allsat: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Cross G++ Linker'
	g++  -o "modsat_allsat" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(C++_DEPS)$(OBJS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) modsat_allsat
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
