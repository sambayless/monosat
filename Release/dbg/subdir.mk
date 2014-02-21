################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../dbg/Debug_mini.cc 

OBJS += \
./dbg/Debug_mini.o 

CC_DEPS += \
./dbg/Debug_mini.d 


# Each subdirectory must supply rules for building sources it contributes
dbg/%.o: ../dbg/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -DNDEBUG -I"/home/sam/workspaceC/modsat/modsat" -I"/home/sam/workspaceC/modsat" -O3 -g3 -Wall -c -fmessage-length=0 -msse2  -std=c++11 -Wno-unused-variable -Wno-unused-but-set-variable -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


