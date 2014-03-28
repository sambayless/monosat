################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Aiger.cpp 

OBJS += \
./Aiger.o 

CPP_DEPS += \
./Aiger.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -DDEBUG_GRAPH -DDEBUG_DIJKSTRA -DDEBUG_MAXFLOW -I/home/sam/workspaceC/modsat/ -I/home/sam/workspaceC/modsat -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11  -Wno-unused-variable -Wno-unused-but-set-variable -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


