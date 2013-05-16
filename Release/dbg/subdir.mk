################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../dbg/Debug.cc \
../dbg/Debug_mini.cc 

OBJS += \
./dbg/Debug.o \
./dbg/Debug_mini.o 

CC_DEPS += \
./dbg/Debug.d \
./dbg/Debug_mini.d 


# Each subdirectory must supply rules for building sources it contributes
dbg/%.o: ../dbg/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -I"/home/sam/workspaceC/modsat" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


