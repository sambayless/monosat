################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../fsm/FSMAcceptDetector.cpp \
../fsm/FSMGeneratesDetector.cpp \
../fsm/FSMGeneratorAcceptorDetector.cpp \
../fsm/FSMTransducesDetector.cpp \
../fsm/P0LAcceptDetector.cpp 

OBJS += \
./fsm/FSMAcceptDetector.o \
./fsm/FSMGeneratesDetector.o \
./fsm/FSMGeneratorAcceptorDetector.o \
./fsm/FSMTransducesDetector.o \
./fsm/P0LAcceptDetector.o 

CPP_DEPS += \
./fsm/FSMAcceptDetector.d \
./fsm/FSMGeneratesDetector.d \
./fsm/FSMGeneratorAcceptorDetector.d \
./fsm/FSMTransducesDetector.d \
./fsm/P0LAcceptDetector.d 


# Each subdirectory must supply rules for building sources it contributes
fsm/%.o: ../fsm/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++1y -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -I.././ -O0 -g3 -Wall -c -fmessage-length=0  -Wno-unused-variable -Wno-unused-but-set-variable   -Wno-sign-compare -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


