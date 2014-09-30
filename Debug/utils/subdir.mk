################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../utils/Options.cc \
../utils/System.cc 

OBJS += \
./utils/Options.o \
./utils/System.o 

CC_DEPS += \
./utils/Options.d \
./utils/System.d 


# Each subdirectory must supply rules for building sources it contributes
utils/%.o: ../utils/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -DDEBUG_MAXFLOW -D__STDC_FORMAT_MACROS -I/home/sam/workspaceC/modsat -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11  -Wno-unused-variable -Wno-unused-but-set-variable  -fsanitize=address -fno-omit-frame-pointer -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


