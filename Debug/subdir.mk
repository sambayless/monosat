################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../Main.cc 

OBJS += \
./Main.o 

CC_DEPS += \
./Main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -DDEBUG_MAXFLOW -D__STDC_FORMAT_MACROS -I/home/sam/workspaceC/modsat -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11  -Wno-unused-variable -Wno-unused-but-set-variable -DDEBUG_GRAPH -DDEBUG_DIJKSTRA -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


