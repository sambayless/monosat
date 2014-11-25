################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../schedule/EDFDetector.cpp 

OBJS += \
./schedule/EDFDetector.o 

CPP_DEPS += \
./schedule/EDFDetector.d 


# Each subdirectory must supply rules for building sources it contributes
schedule/%.o: ../schedule/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -DNDEBUG -DDEBUG_ANALYSIS -I.././ -O3 -g3 -Wall -c -fmessage-length=0 -msse2  -std=c++11 -Wno-unused-variable -Wno-unused-but-set-variable  -static  -static-libgcc     -static-libstdc++ -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


