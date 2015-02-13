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
	clang++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -I.././ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -Wno-missing-braces  -Wno-unused-variable   -Wno-sign-compare -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


