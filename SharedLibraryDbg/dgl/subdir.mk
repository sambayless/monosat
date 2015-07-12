################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../dgl/AllPairs.cpp \
../dgl/ConnectedComponents.cpp \
../dgl/Distance.cpp \
../dgl/Reach.cpp 

OBJS += \
./dgl/AllPairs.o \
./dgl/ConnectedComponents.o \
./dgl/Distance.o \
./dgl/Reach.o 

CPP_DEPS += \
./dgl/AllPairs.d \
./dgl/ConnectedComponents.d \
./dgl/Distance.d \
./dgl/Reach.d 


# Each subdirectory must supply rules for building sources it contributes
dgl/%.o: ../dgl/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -DNDEBUG -I.././ -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -Wno-unused-variable -Wno-unused-but-set-variable  -static  -static-libgcc     -static-libstdc++ -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


