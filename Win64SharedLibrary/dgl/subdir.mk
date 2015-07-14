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
	x86_64-w64-mingw32-g++ -std=c++0x -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -DNDEBUG -I.././ -I.././ -O3 -g3 -Wall -c -fmessage-length=0 -Wno-unused-variable -Wno-unused-but-set-variable -static  -static-libgcc     -static-libstdc++ -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


