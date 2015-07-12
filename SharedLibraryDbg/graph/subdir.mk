################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../graph/AllPairsDetector.cpp \
../graph/ConnectedComponentsDetector.cpp \
../graph/CycleDetector.cpp \
../graph/DistanceDetector.cpp \
../graph/MSTDetector.cpp \
../graph/MaxflowDetector.cpp \
../graph/ReachDetector.cpp \
../graph/SteinerDetector.cpp 

OBJS += \
./graph/AllPairsDetector.o \
./graph/ConnectedComponentsDetector.o \
./graph/CycleDetector.o \
./graph/DistanceDetector.o \
./graph/MSTDetector.o \
./graph/MaxflowDetector.o \
./graph/ReachDetector.o \
./graph/SteinerDetector.o 

CPP_DEPS += \
./graph/AllPairsDetector.d \
./graph/ConnectedComponentsDetector.d \
./graph/CycleDetector.d \
./graph/DistanceDetector.d \
./graph/MSTDetector.d \
./graph/MaxflowDetector.d \
./graph/ReachDetector.d \
./graph/SteinerDetector.d 


# Each subdirectory must supply rules for building sources it contributes
graph/%.o: ../graph/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -DNDEBUG -I.././ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -Wno-unused-variable -Wno-unused-but-set-variable  -static  -static-libgcc     -static-libstdc++ -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


