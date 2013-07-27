################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Aiger.cpp 

CC_SRCS += \
../Main.cc 

OBJS += \
./Aiger.o \
./Main.o 

CC_DEPS += \
./Main.d 

CPP_DEPS += \
./Aiger.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -DDEBUG_SOLVER -DDEBUG_GRAPH -DDEBUG_MAXFLOW -DDEBUG_DIJKSTRA -I"/home/sam/workspaceC/modsat" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

%.o: ../%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -DDEBUG_SOLVER -DDEBUG_GRAPH -DDEBUG_MAXFLOW -DDEBUG_DIJKSTRA -I"/home/sam/workspaceC/modsat" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


