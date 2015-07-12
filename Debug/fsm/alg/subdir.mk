################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../fsm/alg/test.cpp 

OBJS += \
./fsm/alg/test.o 

CPP_DEPS += \
./fsm/alg/test.d 


# Each subdirectory must supply rules for building sources it contributes
fsm/alg/%.o: ../fsm/alg/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -I.././ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11  -Wno-unused-variable -Wno-unused-but-set-variable   -Wno-sign-compare -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


